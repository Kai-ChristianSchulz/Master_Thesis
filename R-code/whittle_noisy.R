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
      lower = c(0.001,0.001),
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
# Parallelize whittle(), apply it to pre-computed daily 1-Second, 1-Minute and 
# 5-Minute realized volatility data from "noisy_estimates.R". 
################################################################################
################################################################################

# We use the "parallel" library for parallelization on a Linux based server.
# For parallelization on a windows based home PC another library should be used
# c.f. https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
library(parallel)

# file path to save the .cvs table containing optimization results
file_path = "whittle_results_noisy.csv"
#file_path = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\whittle_results_noisy.csv"

# Load pre-computed RV data
sigma_hat_1S = scan("RV_1s.txt", sep=" ")
sigma_hat_1 = scan("RV_1Min.txt", sep=" ")
sigma_hat_5 = scan("RV_5Min.txt", sep=" ")

# Combine them into one list
sigma_hat = list("1S" = sigma_hat_1S, "1Min" = sigma_hat_1, "5Min" = sigma_hat_5 )

# Values for m according to each frequency
m = c(86400,1440,288)
names(m) = c("1S","1Min","5Min")


# Initial Optimization values for H and rho
H_init = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
rho_init = c(0.005,0.01,0.1,1,1.5)


# Number of threads for paralellization
#threads =1 
threads = 60


# Creating a grid over which we can paralellize the whittle function
grid = expand.grid(i = 1:length(rho_init), j = 1:length(H_init), freq = names(m), stringsAsFactors = FALSE)

# Create a data frame as placeholder for the optimization results
placeholder = data.frame(
  id = rep(NA, nrow(grid)),
  freq = rep(NA, nrow(grid)),
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

# Save the placeholder as a csv table which we will continously open and save optimization
# results to during paralellization. Each optimization result has a pre-designated
# row for it.
write.csv(placeholder, file = file_path, row.names = FALSE)


# Version of whittle that is fitted for parallelization
paralell_whittle = function(id) {
  
  # Retrieve the index pair from the grid
  i = grid[id, "i"]
  j = grid[id, "j"]
  freq = grid[id, "freq"]
  
  Y = (sigma_hat[[freq]])^2
  
  # Compute the result using whittle estimation
  result = whittle(m[freq], Y, c(H_init[j], rho_init[i]))
  
  # Creates a data frame row where optimization results are saved in
  result_row = data.frame(
    id = id,
    freq = freq,
    m = unname(m[freq]),
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
  
  # Temporary open CSV table
  csv_temp = read.csv(file_path, stringsAsFactors = FALSE)
  
  # Adding result_row in pre-designated row, indicated by id
  csv_temp[id, ] = result_row
  
  # Write updated data back to CSV table
  write.csv(csv_temp, file = file_path, row.names = FALSE)
}

# Run the parallel computation with pblapply (includes progress bar)
mclapply(1:nrow(grid), paralell_whittle, mc.cores = threads)


################################################################################
################################################################################
# Selecting the final adapted whittle estimator as "minimizer among minimizers"
# and saving them in another .csv table. 
################################################################################
################################################################################

# Load the results from the previous section 
optim_results = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\whittle_results_noisy.csv")

# Omit NA entries
optim_results = na.omit(optim_results)

freque = c("1S","1Min","5Min")

# File path to save the final table containing the adapted whittle estimators
# for 1-Sec, 1-Min and 5-Min realized Volatility
file_path_two = "whittle_estimates_noisy.csv"

# Create a data frame as placeholder for all of the adapted whittle estimators
# for each RV frequency
placeholder = data.frame(
  id = rep(NA, length(freque)),
  freq = rep(NA, length(freque)),
  m = rep(NA, length(freque)),
  H_init = rep(NA, length(freque)),
  rho_init = rep(NA, length(freque)),
  H_opt = rep(NA, length(freque)),
  rho_opt = rep(NA, length(freque)),
  value = rep(NA, length(freque)),
  counts = rep(NA, length(freque)),
  convergence = rep(NA, length(freque)),
  message = rep(NA, length(freque)),
  time_Min = rep(NA, length(freque))
)

write.csv(placeholder, file = file_path_two, row.names = FALSE)

# Iteratively we determine for each RV frequency the final adapted whittle 
# estimator as the minimizer among the different minimizers corresponding to 
# each inital value pair combination.

for(i in 1:length(freque)){
  
  # Retrieve the rows of the minimizers corresponding to the current frequency
  temp = optim_results[optim_results$freq == freque[i],]
  
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

# Read the table with whittle estimates
whittle_table = read.csv(file_path_two)


# Reduced table, also containing nu estimate
whittle_table2 = data.frame(
  freq = whittle_table$freq,
  m = whittle_table$m,
  H = whittle_table$H_opt,
  rho = whittle_table$rho_opt,
  nu = whittle_table$rho_opt /2
)

whittle_table2
