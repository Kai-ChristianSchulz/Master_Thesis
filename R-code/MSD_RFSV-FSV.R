library("hypergeo")
# Kumemrs Function M(a,b;z) = genhypergeo(a,b,z)

# Kummers Function
M = function(a,b,z){
  Re(genhypergeo(a,b,z))
}

# Mean squared displacement for stationary fOU
MSD = function(nu,H,th,D){
  x = nu^2/(th^(2*H)) * gamma(2*H+1) *(1-cosh(th*D)) + nu^2 * D^(2*H)
  y = nu^2 * th * D^(2*H+1) /(2*(2*H+1) ) 
  z = exp(-th*D) * M(2*H+1,2*H+2,th*D) - exp(th*D) *M(2*H+1,2*H+2,-th*D)
  return(x-y*z)
}


# Effect of varying theta
nu =0.343
H = 0.5
Delta = 1:20

colors = c("red","blue","green","yellow2","deeppink")
theta = c(0.2,0.3,0.4,0.5, 0.6)

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\MSD_theta.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(Delta), log(MSD(nu,H,0.1,Delta)), type="l", ylab=bquote(log(MSD)), 
     xlab=expression(log(Delta)), lwd=2)
for ( i in 1:length(colors) ){
  lines(log(Delta), log(MSD(nu,H,theta[i],Delta)), col=colors[i], lwd=2 )
}
legend("topleft", col = c("black","red","blue","green","yellow2","deeppink"),
       pch=16, legend=c(bquote(theta == 0.1),bquote(theta == 0.2),
                        bquote(theta == 0.3),bquote(theta == 0.4),
                        bquote(theta == 0.5),bquote(theta == 0.6)), cex =1.2 )
dev.off()

# Effect of varying nu
th =0.5
H = 0.5
Delta = 1:20

colors = c("red","blue","green","yellow2","deeppink")
nu = c(0.2,0.3,0.4,0.5, 0.6)

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\MSD_nu.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(Delta), log(MSD(0.1,H,th,Delta)), type="l", ylab=bquote(log(MSD)), 
     xlab=expression(log(Delta)), ylim=c(-5,1), lwd=2)
for ( i in 1:length(colors) ){
  lines(log(Delta), log(MSD(nu[i],H,th,Delta)), col=colors[i], lwd=2 )
}
legend("topleft", col = c("black","red","blue","green","yellow2","deeppink"),
       pch=16, legend=c(bquote(nu == 0.1),bquote(nu == 0.2),
                        bquote(nu == 0.3),bquote(nu == 0.4),
                        bquote(nu == 0.5),bquote(nu == 0.6)), cex=0.9 )
dev.off()

# Effect of varying H
th =0.01
nu = 0.343
Delta = 1:20

colors = c("red","blue","green","yellow2","deeppink")
H= c(0.2,0.3,0.4,0.5, 0.6)

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\MSD_H.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(Delta), log(MSD(nu,0.1,th,Delta)), type="l", ylab=bquote(log(MSD)), 
     xlab=expression(log(Delta)),ylim=c(-2.2,0), lwd=2)
for ( i in 1:length(colors) ){
  lines(log(Delta), log(MSD(nu,H[i],th,Delta)), col=colors[i], lwd=2 )
}
legend("topleft", col = c("black","red","blue","green","yellow2","deeppink"),
       pch=16, legend=c(bquote(H== 0.1),bquote(H == 0.2),
                        bquote(H == 0.3),bquote(H == 0.4),
                        bquote(H == 0.5),bquote(H == 0.6)), cex =1.2 )
dev.off()

# Empirical MSD
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
RV_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] )
output_SPX = vector(length= length(Delta))

for (i in 1:length(Delta))
  output_SPX[i] =  m(2,Delta[i],RV_SPX)

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\MSD_RFSV_FSV.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(Delta), log(output_SPX), pch=16, ylab=bquote(log(MSD)),  
     xlab=expression(log(Delta)), main="SPX")
lines(log(Delta), log(MSD(0.343,0.129,0.001,Delta)), col="red", lwd=2)
lines(log(Delta), log(MSD(0.343,0.53,0.5,Delta)), col="blue", lwd=2)
legend("bottomright", col = c("black","red","blue"),
       pch=16, legend=c(bquote(m(2,Delta) ),"RFSV","FSV"), cex =1.2  )
dev.off()
