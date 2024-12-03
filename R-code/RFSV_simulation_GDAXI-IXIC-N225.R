# RFSV_simulation_GDAXI-IXIC-N225.R

################################################################################
# We simulate the path of the volatility under the RFSV model as in 
# "RFSV_simulation.R". We simulate the volatility paths for the indices GDAXI,
# IXIC and N225 and plot the paths of the 5 minute realized volatility for 
# comparison.
################################################################################

# Initializing #################################################################
d = 1/100
n = 5000
N = 1/d * n
theta = 2*10^(-5)
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")


################################################################################
# GDAXI
################################################################################

# Initializing
DAX = OM[OM$Symbol==".GDAXI",]
RV5_DAX = sqrt( DAX$rv5[! DAX$rv5 %in% 0] ) 

# GDAXI parameter estimates
H = 0.098
nu = 0.311
M = mean(log(RV5_DAX))

# Simulating fBM
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.098_N-5e+05.txt", sep=" ")
inc_B = diff(B)

# Simulating RFSV model by Euler-Maruyama scheme
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

# Plotting RFSV simulation
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_GDAXI_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility GDAXI",ylab=expression(sigma[t]), 
     xlab=expression(t), ylim =c(0,0.08))
dev.off()

# Plotting 5 minute realized volatility data
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_GDAXI.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(RV5_DAX, type="l", main="5 Min Realized Volatility GDAXI",
     ylab=expression(hat(sigma[t])), xlab=expression(t), ylim = c(0,0.08) )
dev.off()


################################################################################
# IXIC
################################################################################

# Initializing
NASDAQ = OM[OM$Symbol==".IXIC",]
RV5_NASDAQ = sqrt( NASDAQ$rv5[! NASDAQ$rv5 %in% 0] ) 

# IXIC parameter estimates
H = 0.126
nu = 0.297
M = mean(log(RV5_NASDAQ))

# Simulating fBM
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.126_N-5e+05.txt", sep=" ")
inc_B = diff(B)

# Simulating RFSV model by Euler-Maruyama scheme
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

# Plotting RFSV simulation
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_IXIC_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility IXIC",ylab=expression(sigma[t]), 
     xlab=expression(t))
dev.off()

# Plotting 5 minute realized volatility data
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_IXIC.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(RV5_NASDAQ, type="l", main="5 Min Realized Volatility IXIC",
     ylab=expression(hat(sigma[t])), xlab=expression(t), ylim = c(0,0.08) )
dev.off()


################################################################################
# N225
################################################################################

# Initializing
N225 = OM[OM$Symbol==".N225",]
RV5_N225= sqrt( N225$rv5[! N225$rv5 %in% 0] ) 

# N225 parameter estimates
H = 0.119
nu = 0.308
M = mean(log(RV5_N225))

# Simulating fBM
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.119_N-5e+05.txt", sep=" ")
inc_B = diff(B)

# Simulating RFSV model by Euler-Maruyama scheme
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

# Plotting RFSV simulation
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_N225_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility N225",ylab=expression(sigma[t]), 
     xlab=expression(t), ylim =c(0,0.08))
dev.off()

# Plotting 5 minute realized volatility data
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_N225.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(RV5_N225, type="l", main="5 Min Realized Volatility N225",
     ylab=expression(hat(sigma[t])), xlab=expression(t), ylim = c(0,0.08) )
dev.off()

