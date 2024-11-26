# smoothness-volatility_vs_rv.R

################################################################################
# We try to examine the effect of smoothing by using the realized volatility as 
# proxy for volatility. To this end we simulate the volatility "sigma" and the 
# price process "S". We then calculate the daily 5-minute realized volatility
# "sigma_hat" from the price process. We apply the heuristic linear regression 
# based method for estimating the Hurst parameter H from "smoothness.R" to sigma 
# and sigma_hat respectively. We compare the estimated H values.  
# We repeat this procedure for various small values of H0. 
################################################################################


################################################################################
####### Initializing ###########################################################
################################################################################

# We iterate over chosen small values for H0
H0_vec = c(0.01,0.05,0.1,0.15)

# Paramters chosen to represent adapted Whittle estimates of SPX
nu = 1
M = -5

# Other model paramters
theta = 0.0005
N = 720000 # Total number of simulated points
d = 1/288 # Simulating every 5 minutes
n = N*d # Simulating for N*d = 2500 days

# Container storing smoothnes estimates:
# H0 = true H values used for simulation
# H_sigma = heuristic estimate based on sigma
# H_sigma_hat = heuristic estimate based on sigma_hat
results = data.frame(
  H0 = H0_vec,
  H_sigma = rep(NA,length(H0_vec)),
  H_sigma_hat = rep(NA,length(H0_vec))
)


# Functions and parameters related to heuristic smoothnes estimation
m = function(q, Delta, data) {
  sp = 1 
  m = 0  
  while ( sp <= Delta){
    sigma = data[seq(sp,length(data), by = Delta)] 
    m = m + ( sum( abs ( log(sigma[2:length(sigma)]) - log(sigma[1:(length(sigma)-1)]) )^q ) / (length(sigma)-1)) / Delta
    sp = sp+1
  }
  return(m)
}

q = c(0.5,1,1.5,2,3)
Delta = seq(1:50)


################################################################################
####### Calculation and Plotting ###############################################
################################################################################

# Iterating over each value of H0
for(k in 1:length(H0_vec)){
  
  H = H0_vec[k]
  
  # Pre computed fBM at t=1,2,...,720000 using "fBM_circulant_embedding.R"
  B = scan(paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-",H,"_N-720000.txt") )
  inc_B = diff(B)
  
  # Simulating log_sigma via Euler-Maruyama-Scheme every 5 Minutes
  log_sigma = c(M,numeric(N-1))
  for (i in 1:(N-1)){
    log_sigma[i+1] = log_sigma[i] + d * theta * (M - log_sigma[i]) + d^H * nu * inc_B[i] 
  }
  sigma = exp(log_sigma)
  
  # Plotting
  ymax = max(sigma[seq(1,N,by=288)])
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\sigma_H-",H,".eps"),
             width = 16, height = 9,paper = "special", horizontal = FALSE)
  par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
  plot(sigma[seq(1,N,by=288)], type="l", 
       main = bquote("Simulated Volatility," ~ H[0] ~ "=" ~ .(H)),ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,ymax))
  dev.off()
  
  # Simulating log_S via Euler-Maruyama-Scheme every 5 Minutes
  log_S = c(0,numeric(N-1))
  for(i in 1:(N-1)){
    log_S[i+1] = log_S[i] + sqrt(d) * rnorm(1) * sigma[i]
  }
  S = exp(log_S)
  
  # Estimating daily 5-Minute Realized Variance from log_S
  sigma_hat2 = numeric(n)
  for( i in 1:n){
    log_s = log_S[((i-1)*288+1):(i*288)] # Taking the values of log_S from one day
    sigma_hat2[i] = sum( (diff(log_s))^2 ) 
  }
  sigma_hat = sqrt(sigma_hat2)
  
  # Plotting
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\sigma_hat_H-",H,".eps"),
             width = 16, height = 9,paper = "special", horizontal = FALSE)
  par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
  plot(sigma_hat, type="l", 
       main = bquote("Calculated daily 5-min. realized volatility," ~ H[0] ~ "=" ~ .(H)),ylab=expression(hat(sigma[t])), xlab=expression(t), ylim = c(0,ymax))
  dev.off()
  
  # Estimating H for sigma
  # We deliberatly chosse every 288th value of sigma (i.e. once per day)
  # to make it more comparable to sigma_hat, which also consists of only 2500
  # entries.
  output = matrix(nrow = length(q), ncol = length(Delta))
  for (i in 1:length(q)) {
    for (j in 1:length(Delta))
      output[i,j] =  m(q[i],Delta[j],sigma[seq(1,N,by=288)])
  }
  a = vector()
  for (i in 1:length(q))
    a = append(a,lm(log(output[i,])  ~ log(Delta))$coefficients[2] )
  
  # Saving estimate 
  results$H_sigma[k] = unname( lm(a ~ q)$coefficients[2] )
  
  
  # Estimating H for sigma_hat
  output = matrix(nrow = length(q), ncol = length(Delta))
  for (i in 1:length(q)) {
    for (j in 1:length(Delta))
      output[i,j] =  m(q[i],Delta[j],sigma_hat)
  }
  a = vector()
  for (i in 1:length(q))
    a = append(a,lm(log(output[i,])  ~ log(Delta))$coefficients[2] )
  
  # Saving estimate
  results$H_sigma_hat[k] = unname( lm(a ~ q)$coefficients[2] )
}

# Display results
round(results, digits=3)
