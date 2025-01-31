# autocovariance_log-log-plot.R

################################################################################
# Creating the log-log Plots of the auto covariance function of
# 1) log (sigma)
# 2) sigma
################################################################################

# Initializing #################################################################
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = c(".SPX",".GDAXI",".IXIC", ".N225")

H = c(0.129,0.098,0.126,0.119)
names(H) = indices

vol = vector(mode="list", length = 4)
names(vol) = indices

log_vol = vector(mode="list", length = 4)
names(log_vol) = indices

for (index in indices){
  data = sqrt ( OM[OM$Symbol==index,]$rv5 )
  data = data[! data %in% 0]
  vol[[index]] = data
  log_vol[[index]] = log(data)
}

Delta_max = 150

Autocov_Vol = matrix(nrow=Delta_max, ncol=length(indices))
colnames(Autocov_Vol) = indices

Autocov_LogVol = matrix(nrow=Delta_max, ncol=length(indices))
colnames(Autocov_LogVol) = indices

# Calculating Autocovariance for each index and all lags
for (index in indices){
  Autocov_Vol[,index] = acf(vol[[index]], type="covariance", lag.max = Delta_max, 
                    plot=FALSE)$acf[2:(Delta_max+1),1,1]
  Autocov_LogVol[,index] = acf(log_vol[[index]], type="covariance", lag.max = Delta_max, 
                               plot=FALSE)$acf[2:(Delta_max+1),1,1]
}


# Creating the log log plots ###################################################

# For Autocov_Vol
mains = c("SPX","GDAXI","IXIC","N225")
names(mains) = indices

for (index in indices){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Autocov_Vol_",mains[index],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
  x = seq(1:Delta_max)
  plot(log(x),log(Autocov_Vol[,index]),type="p", pch=16, xlab=expression(log(Delta)),  
       ylab = expression(log(Cov(sigma[t], sigma[t + Delta]))), main=mains[index])
  abline(lm(log(Autocov_Vol[,index]) ~ log(x)), col="red")
  dev.off()
}

# For Autocov_LogVol
mains = c("SPX","GDAXI","IXIC","N225")
names(mains) = indices

for (index in indices){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Autocov_LogVol_",mains[index],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
  x = seq(1:Delta_max)
  plot(log(x),log(Autocov_LogVol[,index]),type="p", pch=16, xlab=expression(log(Delta)),  
       ylab = expression(log(Cov(log(sigma[t]), log(sigma[t + Delta])))), main=mains[index])
  abline(lm(log(Autocov_LogVol[,index]) ~ log(x)), col="red")
  dev.off()
}

