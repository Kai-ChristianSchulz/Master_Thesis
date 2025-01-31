# smoothness_all_variances.R

################################################################################
# We re-create the log-log plots of m(q.Delta) vs Delta from "smoothness.R" for 
# SPX. In contrast to before, we will not only be using 5 minute realized 
# volatility as our volatility proxy of choice, but repeat for all available
# volatility estimates for SPX in the OM data set.
################################################################################

# Initializing #################################################################

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

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
SPX = OM[OM$Symbol==".SPX",]


# Realized Variance ############################################################

# 5-Minute 
RV5_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] ) 
# Sub Sampled
RV5SS_SPX = sqrt( SPX$rv5_ss[! SPX$rv5_ss %in% 0] ) 
# 10-Minute
RV10_SPX = sqrt( SPX$rv10[! SPX$rv10 %in% 0] ) 
# Sub Sampled
RV10SS_SPX = sqrt( SPX$rv10_ss[! SPX$rv10_ss %in% 0] ) 


# Realized Kernel Variance #####################################################

#Two-Scale / Bartlett Kernel
RK_twoscale_SPX = sqrt( SPX$rk_twoscale[! SPX$rk_twoscale %in% 0] ) 
# Tukey-Hanning(2) Kernel
RK_th2_SPX = sqrt( SPX$rk_th2[! SPX$rk_th2 %in% 0] ) 
# Parzen Kernel
RKParzen_SPX = sqrt( SPX$rk_parzen[! SPX$rk_parzen %in% 0] )

# Median realized Variance (5-Min)
MedRV_SPX = sqrt( SPX$medrv[! SPX$medrv %in% 0] ) 

# Bipower Variation (robus to jumps cf Barndorff) 
BV_SPX = sqrt( SPX$bv[! SPX$bv %in% 0] ) 
# Sub Sampled
BVSS_SPX = sqrt( SPX$bv_ss[! SPX$bv_ss %in% 0] ) 

################################################################################

vol_SPX = list(RV5_SPX,RV5SS_SPX,RV10_SPX,RV10SS_SPX,RK_twoscale_SPX,RK_th2_SPX,
           RKParzen_SPX,MedRV_SPX,BV_SPX, BVSS_SPX)

names(vol_SPX) = c("RV5","RV5SS","RV10","RV10SS","RK_twoscale","RK_th2",
                   "RKParzen", "MedRV", "BV", "BVSS") 

q = c(0.5,1,1.5,2,3)
delta = seq(1:50)
output = list()

# Calculating m(q,Delta) for each index
for(k in 1:length(vol_SPX)){
  dummy = matrix(nrow = length(q), ncol = length(delta))
  for (i in 1:length(q)) {
    for (j in 1:length(delta))
      dummy[i,j] =  m(q[i],delta[j],vol_SPX[[k]])
  }
  output = append(output,list(dummy))
}

# Container for H, nu and y-intersec estimates
H = numeric(length(vol_SPX))
intersec = numeric(length(vol_SPX))
nu = numeric(length(vol_SPX))


# Creating plots and calculating estimates
for(i in 1:length(vol_SPX)){
  
  # Creating log-log plot of m(q,Delta) vs Delta and adding linear regression
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Smoothness_SPX_",names(vol_SPX)[i],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
  plot(log(delta),log(output[[i]][1,]), col="red", type="p", pch=16, xlim=c(0,4), 
       ylim= c(-3,-0.5), ylab=expression(log(m,Delta)),
       xlab=expression(log(Delta)), main=names(vol_SPX)[i])
  abline(lm(log(output[[i]][1,])  ~ log(delta)), col="red")
  points(log(delta),log(output[[i]][2,]), col="blue", pch=16)
  abline(lm(log(output[[i]][2,])  ~ log(delta)), col="blue")
  points(log(delta),log(output[[i]][3,]), col="green", pch=16)
  abline(lm(log(output[[i]][3,])  ~ log(delta)), col="green")
  points(log(delta),log(output[[i]][4,]), col="yellow2", pch=16)
  abline(lm(log(output[[i]][4,])  ~ log(delta)), col="yellow2")
  points(log(delta),log(output[[i]][5,]), col="deeppink1", pch=16)
  abline(lm(log(output[[i]][5,])  ~ log(delta)), col="deeppink1")
  legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
         pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"))
  dev.off()
  
  # Estimating H and y-intersec
  xi = vector()
  for (k in 1:length(q))
    xi = append(xi,lm(log(output[[i]][k,])  ~ log(delta))$coefficients[2] )
  H[i] = lm(xi ~ q)$coefficients[2]
  intersec[i] = lm(xi ~ q)$coefficients[1]
  
  # Estimating nu
  nu[i] = sd( diff(log(vol_SPX[[i]])))
}

# Display rounded estimates
round(H, digits = 3)
round(nu, digits = 3)
round(intersec, digits =3)
