# Simulating RFSV Modell for DAX, NASDAQ 100, Nikkei225

library(fractionalBM)
set.seed(42)

d = 1/10
theta = 0.0005

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")

################################################################################
# DAX
################################################################################
DAX = OM[OM$Symbol==".GDAXI",]
RV5_DAX = sqrt( DAX$rv5[! DAX$rv5 %in% 0] ) 

n = length(RV5_DAX)
N = 1/d * n

H = 0.098
nu = 0.311
M = mean(log(RV5_DAX))

B = c(0,FBM(N,H,t=n))

X = double(length = (N+1))
X[1] = M

for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_GDAXI_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility GDAXI",ylab=expression(sigma[t]), 
     xlab=expression(t), ylim =c(0,0.08))
dev.off()

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_GDAXI.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(seq(0,n-1), RV5_DAX, type="l", main="5 Min Realized Volatility GDAXI",
     ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.08) )
dev.off()


################################################################################
# NASDAQ
################################################################################
NASDAQ = OM[OM$Symbol==".IXIC",]
RV5_NASDAQ = sqrt( NASDAQ$rv5[! NASDAQ$rv5 %in% 0] ) 

n = length(RV5_NASDAQ)
N = 1/d * n

H = 0.126
nu = 0.297
M = mean(log(RV5_NASDAQ))

B = c(0,FBM(N,H,t=n))

X = double(length = (N+1))
X[1] = M

for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_IXIC_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility IXIC",ylab=expression(sigma[t]), 
     xlab=expression(t), ylim =c(0,0.08))
dev.off()

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_IXIC.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(seq(0,n-1), RV5_NASDAQ, type="l", main="5 Min Realized Volatility IXIC",
     ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.08) )
dev.off()

################################################################################
# N225
################################################################################
N225 = OM[OM$Symbol==".N225",]
RV5_N225= sqrt( N225$rv5[! N225$rv5 %in% 0] ) 

n = length(RV5_N225)
N = 1/d * n

H = 0.119
nu = 0.308
M = mean(log(RV5_N225))

B = c(0,FBM(N,H,t=n))

X = double(length = (N+1))
X[1] = M

for (i in 2:N){
  X[i] = X[i-1] + nu * (B[i] - B[i-1]) + theta * d * (M - X[i-1])
}

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\RFSV_Volatility_N225_Simulated.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(exp(X[seq(1,N,by=1/d)]), type="l", 
     main= "Simulated Volatility N225",ylab=expression(sigma[t]), 
     xlab=expression(t), ylim =c(0,0.08))
dev.off()

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Volatility_N225.eps", width = 16, height = 9,paper = "special", horizontal = FALSE)
par(cex.lab = 1.8 , cex.axis = 1.7, cex.main = 2, mar = c(5, 5, 4, 2) + 0.1)
plot(seq(0,n-1), RV5_N225, type="l", main="5 Min Realized Volatility N225",
     ylab=expression(sigma[t]), xlab=expression(t), ylim = c(0,0.08) )
dev.off()

