# For plotNormalDensity()
library(rcompanion)

# Reading Oxford-Man data set and extracting entries for selected indices
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
SPX = OM[OM$Symbol==".SPX",]
DAX = OM[OM$Symbol==".GDAXI",]
NASDAQ = OM[OM$Symbol==".IXIC",]
N225 = OM[OM$Symbol==".N225",]

# Calculating realized Volatility and omnitting 0 entries
vol_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] )
vol_DAX = sqrt( DAX$rv5[! DAX$rv5 %in% 0] )
vol_NASDAQ = sqrt( NASDAQ$rv5[! NASDAQ$rv5 %in% 0] )
vol_N225 = sqrt( N225$rv5[! N225$rv5 %in% 0] )

# List with the volatilitity data for selected indices
indices = c("SPX","DAX","NASDAQ","N225")
vol = list(vol_SPX, vol_DAX, vol_NASDAQ, vol_N225)
names(vol) = indices

# Container for the increments for each Index for various lags delta
delta = c(1,5,25,125)
increments_SPX = vector(mode = "list", length= length(delta))
increments_DAX = vector(mode = "list", length= length(delta))
increments_NASDAQ = vector(mode = "list", length= length(delta))
increments_N225 = vector(mode = "list", length= length(delta))


# Calculating the increments for various lags

#SPX
for( i in 1:length(delta)) {
  del = delta[i]
  sigma_2 = vol_SPX[seq(del+1,length(vol_SPX), by=1)]
  sigma_1 = vol_SPX[seq(1,length(vol_SPX)-del, by=1)]
  increments_SPX[[i]] = log(sigma_2) - log(sigma_1) 
}

#DAX
for( i in 1:length(delta)) {
  del = delta[i]
  sigma_2 = vol_DAX[seq(del+1,length(vol_DAX), by=1)]
  sigma_1 = vol_DAX[seq(1,length(vol_DAX)-del, by=1)]
  increments_DAX[[i]] = log(sigma_2) - log(sigma_1) 
}

#NASDAQ
for( i in 1:length(delta)) {
  del = delta[i]
  sigma_2 = vol_NASDAQ[seq(del+1,length(vol_NASDAQ), by=1)]
  sigma_1 = vol_NASDAQ[seq(1,length(vol_NASDAQ)-del, by=1)]
  increments_NASDAQ[[i]] = log(sigma_2) - log(sigma_1)
}

#N225
for( i in 1:length(delta)) {
  del = delta[i]
  sigma_2 = vol_N225[seq(del+1,length(vol_N225), by=1)]
  sigma_1 = vol_N225[seq(1,length(vol_N225)-del, by=1)]
  increments_N225[[i]] = log(sigma_2) - log(sigma_1)
}


# Plotting the estimated density for each index for various lags
# Then adding the estimated normal density (solid line) and the re-scaled
# normal density fit (dashed line) for delta = 1

col = c("red","blue","green","yellow2")
xre = seq(-4,4,by=0.01)

# SPX 

H_SPX = 0.129 # Estimated H-value for SPX

for(i in 1:length(delta)){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\dens_SPX",delta[i],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
  # Plotting estimated density and normal density fit
  plotNormalDensity(increments_SPX[[i]], col3=col[i], main="SPX",
                    xlab=bquote(Delta == .(delta[i])), lwd=5 )
  # Evaluating re-scaled normal density and adding it to the plot
  yre = dnorm(xre, delta[i]^H_SPX * mean(increments_SPX[[1]]),
              delta[i]^H_SPX * sd(increments_SPX[[1]]))
  lines(xre,yre,col="red", lty="dashed", lwd=5)
  dev.off()
}

# DAX

H_DAX = 0.098

for(i in 1:length(delta)){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\dens_GDAXI",delta[i],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
  plotNormalDensity(increments_DAX[[i]], col3=col[i], main="GDAXI",
                    xlab=bquote(Delta == .(delta[i])), lwd=5 )
  yre = dnorm(xre, delta[i]^H_DAX * mean(increments_DAX[[1]]),
              delta[i]^H_DAX * sd(increments_DAX[[1]]))
  lines(xre,yre,col="red", lty="dashed", lwd=5)
  dev.off()
}

# NASDAQ

H_NASDAQ = 0.126

for(i in 1:length(delta)){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\dens_IXIC",delta[i],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
  plotNormalDensity(increments_NASDAQ[[i]], col3=col[i], main="IXIC",
                    xlab=bquote(Delta == .(delta[i])), lwd=5 )
  yre = dnorm(xre, delta[i]^H_NASDAQ * mean(increments_NASDAQ[[1]]),
              delta[i]^H_NASDAQ * sd(increments_NASDAQ[[1]]))
  lines(xre,yre,col="red", lty="dashed", lwd=5)
  dev.off()
}

# N225

H_N225 = 0.119

for(i in 1:length(delta)){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\dens_N225_",delta[i],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
  plotNormalDensity(increments_N225[[i]], col3=col[i], main="N225",
                    xlab=bquote(Delta == .(delta[i])), lwd=5)
  yre = dnorm(xre, delta[i]^H_N225 * mean(increments_N225[[1]]),
              delta[i]^H_N225 * sd(increments_N225[[1]]))
  lines(xre,yre,col="red", lty="dashed", lwd=5)
  dev.off()
}


# Calculating Standard Deviation of increments for Delta = 1

sd(increments_SPX[[1]])
sd(increments_DAX[[1]])
sd(increments_NASDAQ[[1]])
sd(increments_N225[[1]])

# Calculating Mean of increments for Delta = 1

mean(increments_SPX[[1]])
mean(increments_DAX[[1]])
mean(increments_NASDAQ[[1]])
mean(increments_N225[[1]])


###############################################################################
# Testing for Normality  
###############################################################################

# We test the increments for normality using the Shapiro-Wilk, Anderson-Darling
# and Jarque-Bera Test

library(nortest) # for Anderson-Darling Test
library(tseries) # for Jarque-Bera Test

# SPX
test_SPX = matrix(nrow=length(delta), ncol = 3)
for (i in 1:length(delta)){
  test_SPX[i,1] = shapiro.test(increments_SPX[[i]][1:5000])$p.value
  test_SPX[i,2] = ad.test(increments_SPX[[i]])$p.value
  test_SPX[i,3] = jarque.bera.test(increments_SPX[[i]])$p.value
}
test_SPX

# DAX
test_DAX = matrix(nrow=length(delta), ncol = 3)
for (i in 1:length(delta)){
  test_DAX[i,1] = shapiro.test(increments_DAX[[i]][1:5000])$p.value
  test_DAX[i,2] = ad.test(increments_DAX[[i]])$p.value
  test_DAX[i,3] = jarque.bera.test(increments_DAX[[i]])$p.value
}
test_DAX

# NASDAQ
test_NASDAQ = matrix(nrow=length(delta), ncol = 3)
for (i in 1:length(delta)){
  test_NASDAQ[i,1] = shapiro.test(increments_NASDAQ[[i]][1:5000])$p.value
  test_NASDAQ[i,2] = ad.test(increments_NASDAQ[[i]])$p.value
  test_NASDAQ[i,3] = jarque.bera.test(increments_NASDAQ[[i]])$p.value
}
test_NASDAQ

# N225
test_N225 = matrix(nrow=length(delta), ncol = 3)
for (i in 1:length(delta)){
  test_N225[i,1] = shapiro.test(increments_N225[[i]][1:5000])$p.value
  test_N225[i,2] = ad.test(increments_N225[[i]])$p.value
  test_N225[i,3] = jarque.bera.test(increments_N225[[i]])$p.value
}
test_N225

