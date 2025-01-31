# empirical_autocovariance.R

################################################################################
# Estimating empirical autokovariance for indices SPX, GDAXI, IXIC, N225 
# as a function of delta^2H 
################################################################################

# Initializing #################################################################
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = c(".SPX",".GDAXI",".IXIC", ".N225")

H = c(0.129,0.098,0.126,0.119)
names(H) = indices

log_vol = vector(mode="list", length = 4)
names(log_vol) = indices

# Extracting log-volatilitiy data
for (index in indices){
  data = sqrt ( OM[OM$Symbol==index,]$rv5 )
  data = data[! data %in% 0]
  log_vol[[index]] = log(data)
}


Delta_max = 50 # Maximum lag size

# Matrix that stores in each column the indices and in each row the related
# covariance for each lag
Cov = matrix(nrow=Delta_max, ncol=length(indices))
colnames(Cov) = indices


# Calculating Autocovariance for each index and all lags #######################
for (index in indices){
    Cov[,index] = acf(log_vol[[index]], type="covariance", lag.max = Delta_max, 
                      plot=FALSE)$acf[2:(Delta_max+1),1,1]
}


# Plotting Cov as a function of Delta^2H and adding linear regression ##########
mains = c("SPX, H = 0.129 ","GDAXI, H = 0.098","IXIC, H = 0.126", "N225, H = 0.119")
names(mains) = indices
file_names = c("SPX","GDAXI","IXIC","N225")
names(file_names) = indices

for (index in indices){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Autocov_log_",file_names[index],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
  x = seq(1:Delta_max) ^(2*H[index])
  plot(x,Cov[,index],type="p", pch=16, xlab=expression(Delta^{2 *H}),  
       ylab = expression(Cov(log(sigma[t]), log(sigma[t + Delta]))), main=mains[index] )
  abline(lm(Cov[,index] ~ x), col="red")
  dev.off()
}
