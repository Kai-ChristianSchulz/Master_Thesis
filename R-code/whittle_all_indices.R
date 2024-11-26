# whittle_all_indices.R

################################################################################
# We apply the adapted Whittle estimator to the 5-minute realized volatility
# data of all indices in the Oxford-Man data set. 
################################################################################


################################################################################
################################################################################
# Implementing the adapted Whittle estimator and related functions
################################################################################
################################################################################

# Calculating autocovariance using FFT:
# For a given data set of length n, this returns a vector consisting of 
# the autocovariance function with lags 0,...,n-1 using fast fourier transform.
# As in the Fukasawa paper, our data are assumed to be drawn from a centered
# random variable. We do NOT demean it.
autocovariance_fft = function(data, demean = FALSE) {
  n = length(data)
  # Demean the data 
  if(demean == TRUE){
    data = data - mean(data)
  }
  # Pad the data with zeros to make it periodic
  data_padded = c(data, rep(0, n))
  # Apply fft
  fft_data = fft(data_padded)
  # Calculate power spectrum (multiply by complex conjugate)
  power_spectrum = fft_data * Conj(fft_data)
  # Normalizing by 1/(2 * n^2) because of how fft is implemented in R
  autocov = Re(fft(power_spectrum, inverse = TRUE)) / (2 * n^2)
  # Return only the first n autocovariance values
  return(autocov[1:n])
}

# C_H
C_H = function(H){
  c = gamma(2*H+1)*sin(pi*H)/(2*pi)
  return(c)
}

# a_{H,rho}(tau,psy)
a = function(H,rho,tau,m,psi=10^(-5),J=20){
  summand = function(j){
    c = (-1)^j * tau^(2*j) / (  factorial(2*j) * rho^2 * C_H(H) )
    c2 = psi^(2*j +2*H) / (2*j + 2*H)
    c3 =  psi^(1+2*j+4*H) / (  rho^2 * C_H(H)* m* pi * (1+2*j+4*H) )
    return( c* (c2 - c3) )
  }
  return(  1/(2*pi) * sum (sapply(0:J, summand)) )
}

# A1_{H,rho}(psi)
A1 = function(H,rho,m,psi=10^(-5)){
  c = psi * log(rho^2 *C_H(H)) 
  c2 = psi * (log(psi) -1)* (1-2*H)
  c3 = psi^(2+2*H) / (rho^2 * C_H(H)*m*pi *(2+2*H) )
  return(1/(2*pi) * (c+c2+c3))
}

# A2_{H,rho}(psi)
A2 = function(H,rho,m,y,psi=10^(-5)){
  n = length(y)
  gamma_hat = autocovariance_fft(y)
  c0 = a(H,rho,0,m,psi) *gamma_hat[1] 
  a = sapply(1:(n-1), function(taui) a(H,rho,taui,m,psi) )
  sum = 2* sum(a*gamma_hat[2:n])
  return(1/(2*pi) *(c0+sum) )
}


# d1_{H}(x,lam)
d1 = function(H,x,la){
  return( (2*pi*x+la)^(-3-2*H) + (2*pi*x-la)^(-3-2*H) )
}

# d2_{H}(x,lam)
d2 = function(H,x,la){
  return( 1/(2*pi*(2+2*H)) *( (2*pi*x+la)^(-2-2*H) + (2*pi*x-la)^(-2-2*H) ) )
}

# l(la) 
l = function(la){
  c = 1/pi * (1- cos(la))
  return( c)
}

# g_{H,rho}(la) (approximation)
g = function(H,rho,la,m,K=500){
  c = rho^2 *C_H(H) *(2*(1-cos(la)))^2 
  d1 = sum ( sapply(1:K, function(ki) d1(H,ki,la) ) )
  d2 = 1/2 * (d2(H,K,la) + d2(H,(K+1),la))
  return(c* ( abs(la)^(-3-2*H) +d1+d2  ) +2/m * l(la) )
}

# I_{n}(la,y) - Periodogramm 
I_n = function(la,y){
  n = length(y)
  j = (1:n)
  sum_result = sum (y * exp(1i*j*la) )
  result = (1/(2*pi*n)) * (abs(sum_result))^2
  return(result)
}



# U_{n}(H,rho) 
U = function(H,rho,m,y,psi=10^(-5)){
  # Dummy functions to vectorize I and g. Needed for integrate().
  I_n_wrapper = function(la){
    val = sapply(la,I_n,y=y)
    return(val)
  }
  g_wrapper = function(la){
    val = sapply(la,g, H=H, rho=rho,m=m)
    return(val)
  }
  # Splitting the integrands. Since log-Term is not as strongly oscillating as
  # Term involving the Periodogramm, we do not need as many subdivions. 
  int1 = integrate(function(la) ( log(g_wrapper(la))  ), 
                   psi, pi,subdivisions = 100)$value
  # For the Term involving the Periodogram, we need at least subdivisions = 1000, 
  # otherwise integral approximation error is too large.
  int2 = integrate(function(la) (  I_n_wrapper(la) / g_wrapper(la) ), 
                   psi, pi, subdivisions = 1000)$value
  A1 = A1(H,rho,m,psi)
  A2 = A2(H,rho,m,y,psi)
  return(1/(2*pi) * (int1+int2) +A1+A2)
}

# Adapted Whittle estimation of H and rho.
# y = vector of differences of log-realized variance
# m = number of daily sampling points used to compute realized variance 
# init = initial values for estimating H and rho
# Output: List containing the optimization results from optim() and time to 
# compute the estimator
whittle = function(m,y,init,psi=10^(-5)){
  # Starting system timer
  timing = system.time({
    # Optimize the function U using optim() and box constraints via the 
    # method L-BFGS-B.
    op =optim(
      par = init,
      fn = function(v) U(v[1],v[2],m,y,psi),
      method = "L-BFGS-B",
      lower = c(0.001,0.1),
      upper = c(0.99,5)
    )
  })
  # Add elapsed time (in minutes) to the list of optimization results
  op = c(op, list(time = timing["elapsed"] / 60))
  op = c(list(initial_val = init), op)
  return(op)
}




################################################################################
################################################################################
# Parallelize whittle(), apply it to Oxford-Man data set saving intermediate 
# results.
################################################################################
# We parallelize the whittle() function from the previous section and apply it 
# to the 5-minute realized variance data of each index in the Oxford-Man data
# set. As described in the master thesis, for each index we compute a total of 
# 44 minimizers, corresponding to the various intial value pairs of H and rho
# used in the optimization. We save all of those 44 optimization results for
# each of the 31 indices in a large .csv table with a total of 1364 rows 
# as intermediate results. 
################################################################################
################################################################################

# We use the "parallel" library for parallelization on a Linux based server.
# For parallelization on a windows based home PC another library should be used
# c.f. https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
library(parallel)

# File path to save the .csv table with intermediate results
file_path = "whittle_results_all.csv"

# Getting Oxford-Man data
OM = read.csv("oxfordmanrealizedvolatilityindices.csv")
indices = unique(OM$Symbol)

# Each index has a unique value for m (see master thesis)
m = c(102,72,102,75,90,84,78,102,102,102,102,78,78,102,78,78,72,78,72,75,96,
      102,102,88,78,78,78,66,102,96,102)
names(m) = indices

# Initial Optimization values for H and rho
H_init = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
rho_init = c(0.5,1.5,2.5,3.5)

# Number of threads for paralellization (adjust to your own needs).
threads = 60

# Creating a grid over which we can parallelize the Whittle function.
# It has 44*31 = 1364 entries, corresponding to each of the 44 initial value pairs 
# (H_init,rho_init) for all 31 Indices in the OM data set. Each row of the grid
# correspond to a unique combination of (rho_init,H_init,index). For example
# row number 46 corresponds to rho_init = 1.5, H_init = 0.01 and the second 
# index, namely ".AORD".
grid = expand.grid(i = 1:length(rho_init), j = 1:length(H_init), index = indices, stringsAsFactors = FALSE)

# Create a data frame as placeholder for all of the 1364 optimization results.
# We can identify each (rho_init,H_init,index) combination by the "id"
# corresponding to its row number in the grid
placeholder = data.frame(
  id = rep(NA, nrow(grid)),
  index = rep(NA, nrow(grid)),
  m = rep(NA, nrow(grid)),
  H_init = rep(NA, nrow(grid)),
  rho_init = rep(NA, nrow(grid)),
  H_opt = rep(NA, nrow(grid)),
  rho_opt = rep(NA, nrow(grid)),
  value = rep(NA, nrow(grid)),
  counts = rep(NA, nrow(grid)),
  convergence = rep(NA, nrow(grid)),
  message = rep(NA, nrow(grid)),
  time_Min = rep(NA, nrow(grid))
)

# Save the placeholder as a .csv table which we will continously open and save 
# optimization results to during paralellization. Each optimization result has 
# a pre-designated row for it, corresponding to its "id".
write.csv(placeholder, file = file_path, row.names = FALSE)


# Version of whittle() function that is fitted for parallelization over the
# different ids or grid rows respectively.
paralell_whittle = function(id) {
  
  # Retrieve the index pair from the grid
  i = grid[id, "i"]
  j = grid[id, "j"]
  index = grid[id, "index"]
  
  # Collect the 5-minute log-realized variance increments of the index
  # from the OM data set. We examine only the last 2500 entries.
  Y = tail(diff(log(OM$rv5[OM$Symbol == index & OM$rv5 != 0])), 2500)
  
  # Compute the adapted whittle estimator
  result = whittle(m[index], Y, c(H_init[j], rho_init[i]))
  
  # Create a data frame row where optimization results are saved in
  result_row = data.frame(
    id = id,
    index = index,
    m = unname(m[index]),
    H_init = H_init[j],
    rho_init = rho_init[i],
    H_opt = result$par[1],
    rho_opt = result$par[2],
    value = result$value,
    counts = paste(unname(result$counts), collapse = ";"),
    convergence = result$convergence,
    message = result$message,
    time_Min = unname(result$time)
  )
  
  # Temporary open .csv table
  csv_temp = read.csv(file_path, stringsAsFactors = FALSE)
  
  # Adding result_row in pre-designated row indicated by id
  csv_temp[id, ] = result_row
  
  # Write updated data back to .csv table
  write.csv(csv_temp, file = file_path, row.names = FALSE)
}

# Run the parallel computation with mclapply
mclapply(1:nrow(grid), paralell_whittle, mc.cores = threads)




################################################################################
################################################################################
# Selecting the final adapted whittle estimator as "minimizer among minimizers"
# for each index and saving them in another .csv table. Also we add 
################################################################################
################################################################################

# Load the results from the previous section 
optim_results = read.csv("whittle_results_all.csv")
indices = unique(optim_results$index)

# File path to save the final table containing the adapted whittle estimators
# (and other optimization results) for each index
file_path_two = "whittle_estimates_all.csv"

# Create a data frame as placeholder for all of the adapted whittle estimators
# for each of the 31 indices.
placeholder = data.frame(
  id = rep(NA, length(indices)),
  index = rep(NA, length(indices)),
  m = rep(NA, length(indices)),
  H_init = rep(NA, length(indices)),
  rho_init = rep(NA, length(indices)),
  H_opt = rep(NA, length(indices)),
  rho_opt = rep(NA, length(indices)),
  value = rep(NA, length(indices)),
  counts = rep(NA, length(indices)),
  convergence = rep(NA, length(indices)),
  message = rep(NA, length(indices)),
  time_Min = rep(NA, length(indices))
)

write.csv(placeholder, file = file_path_two, row.names = FALSE)

# Iteratively we determine for each index the final adapted whittle estimator
# as the minimizer among the 44 minimizers corresponding to each inital value
# pair combination.

for(i in 1:length(indices)){
  
  # Retrieve the rows for the 44 minimizers corresponding to the current index
  temp = optim_results[optim_results$index == indices[i],][1:44,]
 
  # Find out which among those creates the smallest value of the objective function
  # i.e. is the minimizer among minimizers.
  min_ind = which.min( unlist(temp$value) )
  
  # Temporary open .csv table
  csv_temp = read.csv(file_path_two, stringsAsFactors = FALSE)
  
  # Add the row corresponding to the minimizer among minimizer to the table
  csv_temp[i, ] = temp[min_ind,]
  
  # Write updated data back to .csv table
  write.csv(csv_temp, file = file_path_two, row.names = FALSE)
}

# Additionally it is interesting to determine the minimum and maximum values for 
# H and rho among all indices
whittle_table = read.csv(file_path_two)
min(whittle_table$H_opt)
max(whittle_table$H_opt)

min(whittle_table$rho_opt)
max(whittle_table$rho_opt)

# For comparison to the Fukasawa paper we add a column for eta estimates as well
# and for comparison to our previous estimates from "smoothness.R" we add 
# a column with nu estimates.

H = whittle_table$H_opt
rho = whittle_table$rho_opt
eta = (1/250)^(-H) * rho
nu = eta/2 *(1/250)^H

whittle_table2 = data.frame(
  index = whittle_table$index,
  m = whittle_table$m,
  H = whittle_table$H_opt,
  rho = whittle_table$rho_opt,
  eta = (1/250)^(-whittle_table$H_opt) * whittle_table$rho_opt,
  nu = whittle_table$rho_opt /2
)

write.csv(whittle_table2, "whittle_H-rho-eta-nu.csv", row.names = FALSE)

