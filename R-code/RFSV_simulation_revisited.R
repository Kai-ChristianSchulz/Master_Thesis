# RFSV_simulation_revisited.R

################################################################################
# Based on our adapted Whittle estimates of H and nu we simulate the volatility
# process "sigma" using the Euler-Maruyama scheme. Afterwards we simulate the 
# log-prices process "log_S" based on sigma, again by Euler-Maruyama scheme. We
# then compute the daily 5 minute realized volatility from our simulated log_S. 
# Finally, we compare the simulated volatility, the simulation based 5 minute
# realized volatility and the 5 minute realized volatility data from the Oxford-
# Man data set. We repeat this for the four inides SPX, GDAXI, IXIC and N225.
################################################################################


################################################################################
####### Initializing ###########################################################
################################################################################

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")

d = 1/288 # Simulating every 5 minutes
n = 2500 # Simulating 2500 
N = 1/d * n # Total number of simulated points

# Reversion speed chosen to be 1/10 of 1/2500.
theta = 4*10^(-5)

# Adapted Whittle estimates for (H,nu) 
H_vec = c(0.049,0.017,0.028,0.016)
nu_vec = c(0.951,1.352,1.064,1.443)
index_vec = c("SPX","GDAXI","IXIC","N225")


################################################################################
####### Simulating and Plotting ################################################
################################################################################

for(k in 1:length(H_vec)){
  
  H = H_vec[k]
  nu = nu_vec[k]
  index = index_vec[k]
  
  # Loading daily 5 minute realized volatitliy data
  RV5 = sqrt( OM[OM$Symbol==paste0(".",index),]$rv5[! OM[OM$Symbol==paste0(".",index),]$rv5 %in% 0] )
  M = mean(log(RV5))
  
  # Loading pre-computed fBM
  B = scan(paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-",H,"_N-720000.txt") )
  inc_B = diff(B)
  
  # Simulating sigma
  log_sigma = c(M,numeric(N-1))
  for (i in 1:(N-1)){
    log_sigma[i+1] = log_sigma[i] + d * theta * (M - log_sigma[i]) + d^H * nu * inc_B[i] 
  }
  sigma = exp(log_sigma)
  
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
  
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_Vola_",index, ".eps"), width = 16, height = 9,paper = "special", horizontal = FALSE)
  par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
  plot(sigma[seq(1,N,by=288)], type="l", 
       main = paste0("Simulated Volatility ",index),ylab=expression(sigma[t]),
       xlab=expression(t), ylim=c(0,0.09))
  dev.off()
  
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_",index, ".eps"), width = 16, height = 9,paper = "special", horizontal = FALSE)
  par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
  plot(sigma_hat, type="l", 
       main = paste0("Simulated daily 5-min. RV ",index),ylab=expression(hat(sigma[t])),
       xlab=expression(t), ylim=c(0,0.09))
  dev.off()
  
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Data_RV_",index, ".eps"), width = 16, height = 9,paper = "special", horizontal = FALSE)
  par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
  plot(tail(RV5,n), type="l",
       main = paste0("Data on daily 5-min. RV ",index),ylab=expression(hat(sigma[t])),
       xlab=expression(t), ylim=c(0,0.09))
  dev.off()
  
}
