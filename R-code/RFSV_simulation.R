# RFSV_simulation.R

################################################################################
# We simulate the path of the volatility under the RFSV model. To this end, we
# use parameter estimations of SPX. Afterwards we compare the simulated path to
# that of the empirical 5 minute realized volatility of SPX. In addition to that
# we briefly examine if the lag between simulated points D = 1/100 is 
# sufficiently small by investigating if the results compared to D = 1/10 
# differ.
################################################################################


################################################################################
# Simulating volatility under the RFSV Model for SPX parameter estimations 
################################################################################

# We set D = 1/100, meaning we simulate 100 points per day. Furthermore we 
# use the SPX parameter estimations H = 0.129, nu = 0.343 and set
# theta = 2*10^(-5), M = X_0 = -4.956, n = 5000.

d = 1/100 # Simulating 100 points per day
n = 5000 # Simulating a total of 5000 days
N = 1/d * n # Total number of points to be simulated (half a million)

# Estimated parameters
H = 0.129 
nu = 0.343
theta = 2*10^(-5)
M = -4.956 # Estimated as mean of 5-Min Realized log-Volatility

# Since Simulation is quite long, we have pre-computed fBM values
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.129_N-5e+05.txt", sep=" ")
inc_B = diff(B)

# Simulate RFSV model by Euler-Maruyama scheme
# Note that B[1], B[2] are generated at t=1,2,... Therefore we need to re-scale
# with d^H.
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

# Plotting RFSV simulation
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# Plotting 5 minute realized volatility data of SPX from OM data set
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
SPX = OM[OM$Symbol==".SPX",]
RV5_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] ) 

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_SPX.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(RV5_SPX, type="l", main="5 Min Realized Volatility SPX",
     ylab=expression(hat(sigma[t])), xlab=expression(t), ylim = c(0,0.09) )
dev.off()



################################################################################
# Simulating RFSV for D = 1/10 and D = 1/100
################################################################################

# D = 1/10 #####################################################################

d = 1/10
N = 1/d * n

# First simulation
X = double(length = N)
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.129_N-50000.txt", sep=" ")
inc_B = diff(B)
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated2.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# Second simulation
X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.129_N-50000_2.txt", sep=" ")
inc_B = diff(B)
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated3.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# D = 1/100 ####################################################################

d = 1/100
N = 1/d * n

# First simulation
X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.129_N-5e+05_2.txt", sep=" ")
inc_B = diff(B)
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated4.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# Second simulation
X = double(length = (N+1))
inc_B = diff(B)
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-0.129_N-5e+05_3.txt", sep=" ")
X = c(M,numeric(N-1))
for (i in 1:(N-1)){
  X[i+1] = X[i] + nu *d^H* inc_B[i] + theta * d * (M - X[i])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated5.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()





