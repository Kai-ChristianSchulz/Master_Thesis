# We investigate the influence of the approximation error of using realized 
# instead of spot variance on our estimation of the Hurst paramter via linear
# regression when the signal is weak.
# We do so by simulating spot variance and afterwards estimating 1Second, 1Min 
# and 5Min realized variance from it.


################################################################################
# 1) Simulating sigma_t 
################################################################################

set.seed(2024)

d = 1/86400 # Simulating every second
n = 2000  # Simulating n=1000 (trading) days total
N = 1/d * n  # Number of total simulated points

nu = 0.01
M = -1.6
theta = 0.001

# Simualting log-volatility X
X = double(length = N )
X[1] = M

for (i in 2:length(X)){
  X[i] = X[i-1] + nu * sqrt(d) * rnorm(1) + theta * d * (M - X[i-1])
}

# Getting volatility
sigma = exp(X)

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Simulated_Spot_Vola.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(sigma[seq(1,length(sigma), by=1/d)], type="l", xlab="trading days",
     ylab="Volatility", main="Simulated (Spot) Volatility" )
dev.off()



################################################################################
# 2) Simulating P_t
################################################################################

# Simulating price path
S = double(length = N)
S[1] = 1400

for ( i in 2:N){
  S[i] = S[i-1] + sigma[i-1] * S[i-1] * sqrt(d) * rnorm(1) 
}

#plot(S[seq(1,length(S), by=1/d)], type="l")

################################################################################
# 3) Compute Realized Variance
################################################################################

# 1-Second RV 
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_1S = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_1S,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_1s.txt")

# 1-Min RV 
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  s = s[seq(1,length(s), by =60)] # Taking the value each 60 seconds
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_1 = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_1,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_1Min.txt")

# 5-Min RV
sigma2_hat = double(length = n )

for ( i in 1:length(sigma2_hat)){
  s = S[((i-1)*1/d+1):(i*1/d)] # Prices for one day
  s = s[seq(1,length(s), by =300)] # Taking the value each 5 Minutes
  sigma2_hat[i] = sum ( (diff(log(s)))^2 )
}

sigma_hat_5 = sqrt(sigma2_hat[! sigma2_hat %in% 0])
write(sigma_hat_5,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_5Min.txt")

################################################################################
# 4) Plot pre-computed Realized Variance
################################################################################

# Load pre-computed RV data
sigma_hat_1S = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_1s.txt", sep=" ")
sigma_hat_1 = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_1Min.txt", sep=" ")
sigma_hat_5 = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Sim_RV_5Min.txt", sep=" ")

# 1S RV
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_1s.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(sigma_hat_1S, type="l", xlab="trading days",
     ylab="Volatility", main="Simulated 1 Second Realized Volatility")
dev.off()

# 1Min RV
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_1Min.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(sigma_hat_1, type="l",xlab="trading days",
     ylab="Volatility", main="Simulated 1 Minute Realized Volatility")
dev.off()

# 5Min RV
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Sim_RV_5Min.eps", width = 6, height = 6,paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(sigma_hat_5, type="l",xlab="trading days",
     ylab="Volatility", main="Simulated 5 Minute Realized Volatility")
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

# Original Volatility
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    outputOG[i,j] =  m(q[i],delta[j],sigma[seq(1,length(sigma), by=1/d)])
}

# 1S
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output1S[i,j] =  m(q[i],delta[j],sigma_hat_1S)
}

# 1 Min 
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output1[i,j] =  m(q[i],delta[j],sigma_hat_1)
}

# 5 Min
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output5[i,j] =  m(q[i],delta[j],sigma_hat_5)
}


# Plot for 5 Min RV
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

# Linear Regression 
aOG = vector()
for (i in 1:length(q))
  aOG = append(aOG,lm(log(outputOG[i,])  ~ log(delta))$coefficients[2] )


a1S = vector()
for (i in 1:length(q))
  a1S = append(a1S,lm(log(output1S[i,])  ~ log(delta))$coefficients[2] )

a1 = vector()
for (i in 1:length(q))
  a1 = append(a1,lm(log(output1[i,])  ~ log(delta))$coefficients[2] )

a5 = vector()
for (i in 1:length(q))
  a5 = append(a5,lm(log(output5[i,])  ~ log(delta))$coefficients[2] )

# Plotting Regression
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

# H estimates
lm(aOG ~ q)$coefficients[2] 
lm(a1S ~ q)$coefficients[2] 
lm(a1 ~ q)$coefficients[2] 
lm(a5 ~ q)$coefficients[2] 