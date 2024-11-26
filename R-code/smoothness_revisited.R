# smoothness_revisited.R

################################################################################
# We simulated volatility "sigma" and price process "S" based on the adapted 
# Whittle estimates for H and nu. Then, we compute the daily 5-minute 
# realized volatility  "sigma_hat" from the simulated price process. We apply 
# the heuristic estimation methods for H and nu from "smoothness.R" to sigma_hat.
################################################################################


################################################################################
####### Initializing ###########################################################
################################################################################

# Pre computed fBM at t=1,2,...,720000
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.049_N-720000.txt")
inc_B = diff(B)

# Parameters
H = 0.049
nu = 0.951

M = -4.956
theta = 0.0005

N = length(B)
d = 1/288 # Corresponding to the length of 5 Minutes in a day of 24 hours
n = N*d

################################################################################
####### Simulating sigma, S and computing sigma_hat ############################
################################################################################

# Simulating log_sigma via Euler-Maruyama-Scheme every 5 Minutes
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

################################################################################
####### Estimating H and nu via heuristics #####################################
################################################################################

# To estimate H we use the same two linear regression approaches as in
# "smoothnes.R"

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

output = matrix(nrow = length(q), ncol = length(Delta))
for (i in 1:length(q)) {
  for (j in 1:length(Delta))
    output[i,j] =  m(q[i],Delta[j],sigma_hat)
}

# Plotting first linear regression
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\smoothness_revisited.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(Delta),log(output[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylim= c(-3,0), ylab=expression(log(m,Delta)),
     xlab=expression(log(Delta)))
abline(lm(log(output[1,])  ~ log(Delta)), col="red")
points(log(Delta),log(output[2,]), col="blue", pch=16)
abline(lm(log(output[2,])  ~ log(Delta)), col="blue")
points(log(Delta),log(output[3,]), col="green", pch=16)
abline(lm(log(output[3,])  ~ log(Delta)), col="green")
points(log(Delta),log(output[4,]), col="yellow2", pch=16)
abline(lm(log(output[4,])  ~ log(Delta)), col="yellow2")
points(log(Delta),log(output[5,]), col="deeppink1", pch=16)
abline(lm(log(output[5,])  ~ log(Delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex= 1.2)
dev.off()


a = vector()
for (i in 1:length(q))
  a = append(a,lm(log(output[i,])  ~ log(Delta))$coefficients[2] )

# Plotting second linear regression
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\a_revisited.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,a, type="p", pch=16, xlim=c(0,3), ylim=c(0,0.5),xlab=expression(q), 
     ylab=expression(a[q]))
abline(lm(a~ q), col="red")
dev.off()

# Estimating H as slope of a(q)
H_est = unname( lm(a ~ q)$coefficients[2] )
H_est

# Estimating nu as the empirical standard deviation of log(sigma_hat) increments
nu_est = sd( diff(log(sigma_hat)))
nu_est
