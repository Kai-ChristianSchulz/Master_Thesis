# noisy_estimates.R

################################################################################
# We simulate volatility having a weak signal nu^2 / (2*theta) and estimate
# daily 1-Sec, 1-Min and 5-Min realized volatility from it. Afterwards we apply
# the same linear-regression-based method as in "smoothness.R" to estimate
# the Hurst parameter H of the volatility and realized volatilities.
################################################################################

################################################################################
# Initializing 
################################################################################

set.seed(42)
d = 1/86400  # Simulating every second
n = 2500  # Simulating n=2500 (trading) days total
N = 1/d * n  # Number of total simulated points

# We intentionally choose a weak signal nu^2 / (2*theta) = 0.05
nu = 0.01
M = -5
theta = 0.001

################################################################################
# 1) Simulating volatility
################################################################################

# Simulating log_volatility via Euler-Maruyama-Scheme 
log_sigma = c(M,numeric(N-1))
for (i in 1:(N-1)){
  log_sigma[i+1] = log_sigma[i] + d * theta * (M - log_sigma[i]) + sqrt(d)* nu * rnorm(1) 
}
sigma = exp(log_sigma)

# Plotting volatility sigma 
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_Vola.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
plot(sigma[seq(1,length(sigma), by=1/d)], type="l", xlab=expression(t),
     ylab=expression(sigma[t]), main="Simulated Volatility",ylim = c(0.005,0.0085) )
dev.off()

################################################################################
# 2) Simulating price process
################################################################################

# Simulating log_S via Euler-Maruyama-Scheme 
log_S = c(0,numeric(N-1))
for(i in 1:(N-1)){
  log_S[i+1] = log_S[i] + sqrt(d) * rnorm(1) * sigma[i]
}
S = exp(log_S)

################################################################################
# 3) Compute Realized Variances
################################################################################

# 1-Second RV ##################################################################
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_1S = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_1S,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_1s.txt")

# 1-Min RV #####################################################################
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  s = s[seq(1,length(s), by =60)] # Taking the value each 60 seconds
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_1 = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_1,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_1Min.txt")

# 5-Min RV #####################################################################
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  s = s[seq(1,length(s), by =300)] # Taking the value each 5 Minutes
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_5 = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_5,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_5Min.txt")


################################################################################
# 4) Plot (pre-computed) realized volatility
################################################################################

# Load pre-computed RV data
sigma_hat_1S = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_1s.txt", sep=" ")
sigma_hat_1 = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_1Min.txt", sep=" ")
sigma_hat_5 = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\RV_5Min.txt", sep=" ")

# 1S RV ########################################################################
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_1s.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
plot(sigma_hat_1S, type="l", ylab=expression(hat(sigma[t])),
     xlab=expression(t), main="Simulated daily 1-sec. RV",ylim = c(0.005,0.0085))
dev.off()

# 1Min RV ######################################################################
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_1Min.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
plot(sigma_hat_1, type="l",ylab=expression(hat(sigma[t])),
     xlab=expression(t), main="Simulated daily 1-min. RV",ylim = c(0.005,0.0085))
dev.off()

# 5Min RV ######################################################################
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_5Min.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
plot(sigma_hat_5, type="l",ylab=expression(hat(sigma[t])),
     xlab=expression(t), main="Simulated daily 5-min. RV",ylim = c(0.005,0.0085))
dev.off()


################################################################################
# 5) Estimating H via linear regression approach
################################################################################

m = function(q, delta, data) {
  sp = 1 # Starting Position
  m = 0  # Returned estimate
  while ( sp <= delta){
    sigma = data[seq(sp,length(data), by = delta)] # choosing grid points from data set
    # Calculating the estimate (already avering)
    m = m + ( sum( abs ( log(sigma[2:length(sigma)]) - log(sigma[1:(length(sigma)-1)]) )^q ) / (length(sigma)-1)) / delta
    sp = sp+1 # Shifting starting position for grid points
  }
  return(m)
}


q = c(0.5,1,1.5,2,3)
delta = seq(1:50)

# Output matrix contains calculations for m(q,delta)
# Rows: different q values
# Cols: different delta values
outputOG = matrix(nrow = length(q), ncol = length(delta))
output1S = matrix(nrow = length(q), ncol = length(delta))
output1 = matrix(nrow = length(q), ncol = length(delta))
output5 = matrix(nrow = length(q), ncol = length(delta))

# Volatility 
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    outputOG[i,j] =  m(q[i],delta[j],sigma[seq(1,length(sigma), by=1/d)])
}

# 1S realized volatlity
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output1S[i,j] =  m(q[i],delta[j],sigma_hat_1S)
}

# 1Min realized volatlity
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output1[i,j] =  m(q[i],delta[j],sigma_hat_1)
}

# 5Min realized volatlity
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output5[i,j] =  m(q[i],delta[j],sigma_hat_5)
}


# Creating log-log plot of m(q,Delta) vs Delta for 5Min realized volatlity
# (We do not compute it for 1Min and 1S realized volatlity)
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Simulated_5Min_RV_Regression.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(delta),log(output5[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylab=expression(log(m,Delta)), ylim = c(-11,0),
     xlab=expression(log(Delta)), main="5 Minute Realized Volatility")
abline(lm(log(output5[1,])  ~ log(delta)), col="red")
points(log(delta),log(output5[2,]), col="blue", pch=16)
abline(lm(log(output5[2,])  ~ log(delta)), col="blue")
points(log(delta),log(output5[3,]), col="green", pch=16)
abline(lm(log(output5[3,])  ~ log(delta)), col="green")
points(log(delta),log(output5[4,]), col="yellow2", pch=16)
abline(lm(log(output5[4,])  ~ log(delta)), col="yellow2")
points(log(delta),log(output5[5,]), col="deeppink1", pch=16)
abline(lm(log(output5[5,])  ~ log(delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex =1.1)
dev.off()

# Computing linear regression for volatility and realized volatlity

# volatility
aOG = vector()
for (i in 1:length(q))
  aOG = append(aOG,lm(log(outputOG[i,])  ~ log(delta))$coefficients[2] )

# 1S realized volatlity
a1S = vector()
for (i in 1:length(q))
  a1S = append(a1S,lm(log(output1S[i,])  ~ log(delta))$coefficients[2] )

# 1Min realized volatlity
a1 = vector()
for (i in 1:length(q))
  a1 = append(a1,lm(log(output1[i,])  ~ log(delta))$coefficients[2] )

# 5Min realized volatlity
a5 = vector()
for (i in 1:length(q))
  a5 = append(a5,lm(log(output5[i,])  ~ log(delta))$coefficients[2] )

# Plotting a_q vs q and adding linear regression for all volatilities
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Simulated_Regression.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,aOG, type="p", pch=16, xlim=c(0,3),xlab=expression(q), 
     ylab=expression(a[q]), ylim= c(0,1.7))
abline(lm(aOG ~ q), col="black", lwd=2)
points(q,a1S, pch= 16)
abline(lm(a1S ~ q), col="red", lwd=2)
points(q,a1, pch= 16)
abline(lm(a1 ~ q), col="blue", lwd=2)
points(q,a5, pch= 16)
abline(lm(a5 ~ q), col="green", lwd=2)
legend("topleft", col = c("black","red","blue","green"), lty=c(1,1,1,1), lwd=2,
       legend=c("Spot","1 Sec","1 Min","5 Min"),cex=1.2)
dev.off()

# Displaying H estimates
lm(aOG ~ q)$coefficients[2] 
lm(a1S ~ q)$coefficients[2] 
lm(a1 ~ q)$coefficients[2] 
lm(a5 ~ q)$coefficients[2] 
