################################################################################
# Simulating Volatility under the RFSV Model
################################################################################

# We put delta = 1/100 meaning we simulate 100 points per day
# H = 0.129, nu = 0.343, theta = 0.0005, M = X_0 = -5
# n = 5052 

d = 1/100
n = 5052
N = 1/d * n

H = 0.129 
nu = 0.343
theta = 0.0005
M = -4.956 # Estimated as mean of 5-Min Realized log-Volatility

# Since Simulation is quite long, we have pre-computed fBM values
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation.txt", sep=" ")

# Computing Modell
X = double(length = (N+1))
X[1] = M

for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

# Plot simulated Model

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# Plot real 5 Minute Realized Volatility Data
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
SPX = OM[OM$Symbol==".SPX",]
RV5_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] ) 

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_SPX.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(seq(0,5051), RV5_SPX, type="l", main="5 Min Realized Volatility SPX",
     ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09) )
dev.off()

################################################################################
# Plotting Volatility for the first 50 days
################################################################################
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_50days.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[1:(50/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), 
     xaxt = 'n', xlab= "t", ylim = c(0,0.09))
axis(side=1, at = c(0,1000,2000,3000,4000,5000) , labels =c(0,10,20,30,40,50))
dev.off()


################################################################################
# Plotting other simulations for various delta
################################################################################

# d=1/10 #######################################################################

d = 1/10
N = 1/d * n
X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation2.txt", sep=" ")
X[1] = M
for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated2.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation3.txt", sep=" ")
X[1] = M
for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated3.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

# d = 1/100 ###################################################################

d = 1/100
N = 1/d * n
X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation4.txt", sep=" ")
X[1] = M
for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated4.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()

X = double(length = (N+1))
B = scan("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation5.txt", sep=" ")
X[1] = M
for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_SPX_Simulated5.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility SPX",ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.09))
dev.off()





