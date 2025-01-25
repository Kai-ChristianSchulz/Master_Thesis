# smoothness.R

################################################################################
# We investigate several properties of the volatility process. To this
# end, we use the 5 minute realized volatility from the OM data set as an 
# empirical proxy. 
################################################################################


################################################################################
# Calculating and plotting m(q,delta)
################################################################################
# We calculate an empirical estimate of the q-th absolute moment of the 
# log increments of the volatility process, i.e. E(|log(signa_delta) - log(sigma_0) |^q)
# We do this for the indices SPX, GDAXI, IXIC and N225 using 5 minute realized
# volatility data from the OM data set.
################################################################################

# Empirical q-th absolute moment of log-increments
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


# Reading Oxford-Man data set and extracting entries for selected indices 
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
SPX = OM[OM$Symbol==".SPX",]
DAX = OM[OM$Symbol==".GDAXI",]
NASDAQ = OM[OM$Symbol==".IXIC",]
N225 = OM[OM$Symbol==".N225",]

# Calculating realized Volatility and omnitting 0 entries
RV_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] )
RV_DAX = sqrt( DAX$rv5[! DAX$rv5 %in% 0] )
RV_NASDAQ = sqrt( NASDAQ$rv5[! NASDAQ$rv5 %in% 0] )
RV_N225 = sqrt( N225$rv5[! N225$rv5 %in% 0] )

# Calculating m(q,delta) for each indice
q = c(0.5,1,1.5,2,3)
delta = seq(1:50)
output_SPX = matrix(nrow = length(q), ncol = length(delta))
output_DAX = matrix(nrow = length(q), ncol = length(delta))
output_NASDAQ = matrix(nrow = length(q), ncol = length(delta))
output_N225 = matrix(nrow = length(q), ncol = length(delta))

# For SPX
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output_SPX[i,j] =  m(q[i],delta[j],RV_SPX)
}

# For GDAXI
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output_DAX[i,j] =  m(q[i],delta[j],RV_DAX)
}

# For IXIC
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output_NASDAQ[i,j] =  m(q[i],delta[j],RV_NASDAQ)
}

# For N225
for (i in 1:length(q)) {
  for (j in 1:length(delta))
    output_N225[i,j] =  m(q[i],delta[j],RV_N225)
}

# Plotting m(q,delta) and fitted linear regression #############################

# For SPX
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Smoothness_SPX.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(delta),log(output_SPX[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylim= c(-3,-0.5), ylab=expression(log(m(q,Delta))),
     xlab=expression(log(Delta)), main="SPX")
abline(lm(log(output_SPX[1,])  ~ log(delta)), col="red")
points(log(delta),log(output_SPX[2,]), col="blue", pch=16)
abline(lm(log(output_SPX[2,])  ~ log(delta)), col="blue")
points(log(delta),log(output_SPX[3,]), col="green", pch=16)
abline(lm(log(output_SPX[3,])  ~ log(delta)), col="green")
points(log(delta),log(output_SPX[4,]), col="yellow2", pch=16)
abline(lm(log(output_SPX[4,])  ~ log(delta)), col="yellow2")
points(log(delta),log(output_SPX[5,]), col="deeppink1", pch=16)
abline(lm(log(output_SPX[5,])  ~ log(delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex= 1.2)
dev.off()

# For GDAXI
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Smoothness_GDAXI.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(delta),log(output_DAX[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylim= c(-3,-0.5), ylab=expression(log(m(q,Delta))),
     xlab=expression(log(Delta)), main="GDAXI")
abline(lm(log(output_DAX[1,])  ~ log(delta)), col="red")
points(log(delta),log(output_DAX[2,]), col="blue", pch=16)
abline(lm(log(output_DAX[2,])  ~ log(delta)), col="blue")
points(log(delta),log(output_DAX[3,]), col="green", pch=16)
abline(lm(log(output_DAX[3,])  ~ log(delta)), col="green")
points(log(delta),log(output_DAX[4,]), col="yellow2", pch=16)
abline(lm(log(output_DAX[4,])  ~ log(delta)), col="yellow2")
points(log(delta),log(output_DAX[5,]), col="deeppink1", pch=16)
abline(lm(log(output_DAX[5,])  ~ log(delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex= 1.2)
dev.off()

# For IXIC
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Smoothness_IXIC.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(delta),log(output_NASDAQ[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylim= c(-3,-0.5), ylab=expression(log(m(q,Delta))),
     xlab=expression(log(Delta)), main="IXIC")
abline(lm(log(output_NASDAQ[1,])  ~ log(delta)), col="red")
points(log(delta),log(output_NASDAQ[2,]), col="blue", pch=16)
abline(lm(log(output_NASDAQ[2,])  ~ log(delta)), col="blue")
points(log(delta),log(output_NASDAQ[3,]), col="green", pch=16)
abline(lm(log(output_NASDAQ[3,])  ~ log(delta)), col="green")
points(log(delta),log(output_NASDAQ[4,]), col="yellow2", pch=16)
abline(lm(log(output_NASDAQ[4,])  ~ log(delta)), col="yellow2")
points(log(delta),log(output_NASDAQ[5,]), col="deeppink1", pch=16)
abline(lm(log(output_NASDAQ[5,])  ~ log(delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex= 1.2)
dev.off()

# For N225
postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\Smoothness_N225.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5, mar = c(5, 5, 4, 2) + 0.1)
plot(log(delta),log(output_N225[1,]), col="red", type="p", pch=16, xlim=c(0,4), 
     ylim= c(-3,-0.5), ylab=expression(log(m(q,Delta))),
     xlab=expression(log(Delta)), main="N225")
abline(lm(log(output_N225[1,])  ~ log(delta)), col="red")
points(log(delta),log(output_N225[2,]), col="blue", pch=16)
abline(lm(log(output_N225[2,])  ~ log(delta)), col="blue")
points(log(delta),log(output_N225[3,]), col="green", pch=16)
abline(lm(log(output_N225[3,])  ~ log(delta)), col="green")
points(log(delta),log(output_N225[4,]), col="yellow2", pch=16)
abline(lm(log(output_N225[4,])  ~ log(delta)), col="yellow2")
points(log(delta),log(output_N225[5,]), col="deeppink1", pch=16)
abline(lm(log(output_N225[5,])  ~ log(delta)), col="deeppink1")

legend("bottomright", col = c("red","blue","green","yellow2","deeppink"),
       pch=16, legend=c("q=0.5","q=1","q=1.5","q=2","q=3"), cex= 1.2)
dev.off()

################################################################################
# log-log plot suggests the linear relation: 
# log(m(q,delta)) = a_q * log(delta) + c_q
# Re-formulating, yields:
# m(q,delta) = delta ^a_q * exp(c_q) = delta^a_q * b_q
################################################################################


################################################################################
# Calculating a_q vs q
###############################################################################
# We investigate a_q for different q by plotting a_q vs q and 
# then applying a linear regression. We again do this for SPX, GDAXI, IXIC and 
# N225.
###############################################################################

# SPX 
a_SPX = vector()
for (i in 1:length(q))
  a_SPX = append(a_SPX,lm(log(output_SPX[i,])  ~ log(delta))$coefficients[2] )

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\a_SPX.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,a_SPX, type="p", pch=16, xlim=c(0,3), ylim=c(0,0.4),xlab=expression(q), 
     ylab=expression(a[q]),main="SPX")
abline(lm(a_SPX ~ q), col="red")
dev.off()

lm(a_SPX ~ q) 


# GDAXI
a_DAX = vector()
for (i in 1:length(q))
  a_DAX = append(a_DAX,lm(log(output_DAX[i,])  ~ log(delta))$coefficients[2] )

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\a_GDAXI.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,a_DAX, type="p", pch=16, xlim=c(0,3), ylim=c(0,0.4),xlab=expression(q), 
     ylab=expression(a[q]),main="GDAXI")
abline(lm(a_DAX~ q), col="red")
dev.off()

lm(a_DAX ~ q) 


# IXIC
a_NASDAQ = vector()
for (i in 1:length(q))
  a_NASDAQ = append(a_NASDAQ,lm(log(output_NASDAQ[i,])  ~ log(delta))$coefficients[2] )

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\a_IXIC.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,a_NASDAQ, type="p", pch=16, xlim=c(0,3), ylim=c(0,0.4),xlab=expression(q), 
     ylab=expression(a[q]),main="IXIC")
abline(lm(a_NASDAQ ~ q), col="red")
dev.off()

lm(a_NASDAQ ~ q) 


# N225
a_N225 = vector()
for (i in 1:length(q))
  a_N225 = append(a_N225,lm(log(output_N225[i,])  ~ log(delta))$coefficients[2] )

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\a_N225.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5)
plot(q,a_N225, type="p", pch=16, xlim=c(0,3), ylim=c(0,0.4),xlab=expression(q), 
     ylab=expression(a[q]),main="N225")
abline(lm(a_N225 ~ q), col="red")
dev.off()

lm(a_N225 ~ q) 

###############################################################################
# Calculating b_q vs q  
###############################################################################
# We investigate b_q for different q. To this end, we calculate 
# nu^q * K_q where nu is the (empirical) standard deviation of the log-increments
# (with lag 1) of the volatility process and K_q is the q-th absolute moment
# of standard Gaussian variable. Again, repeating for SPX, GDAXI, IXIC and 
# N225 respectively.
###############################################################################

# Calculating bq for each q and each Index
bq_SPX = vector(length= length(q))
for (i in 1:length(q))
  bq_SPX[i] = exp(lm(log(output_SPX[i,])  ~ log(delta))$coefficients[1]) 

bq_DAX = vector(length= length(q))
for (i in 1:length(q))
  bq_DAX[i] = exp(lm(log(output_DAX[i,])  ~ log(delta))$coefficients[1]) 

bq_NASDAQ = vector(length= length(q))
for (i in 1:length(q))
  bq_NASDAQ[i] = exp(lm(log(output_NASDAQ[i,])  ~ log(delta))$coefficients[1]) 

bq_N225 = vector(length= length(q))
for (i in 1:length(q))
  bq_N225[i] = exp(lm(log(output_N225[i,])  ~ log(delta))$coefficients[1]) 


# q-th absolute moment of a standard Gaussian variable (for q=0.5,1,1.5,2,3)
Kq = c(0.822,0.798,0.860,1,1.596)

# Empirical Standard Deviation for the increments log sigma_t+1 - log sigma
# for each index. Values calculated in "Distribution Increments.R"
nu_SPX = 0.3426148
nu_DAX = 0.3113694
nu_NASDAQ = 0.2970392
nu_N225 = 0.3078223

# Calculating nu^q * Kq for each q and comparing it to bq
est_SPX = cbind(q,bq_SPX , nu_SPX ^q * Kq, bq_SPX - nu_SPX ^q *Kq) 
est_DAX = cbind(q,bq_DAX, nu_DAX ^q * Kq, bq_DAX - nu_DAX ^q * Kq)
est_NASDAQ = cbind(q,bq_NASDAQ, nu_NASDAQ ^q * Kq, bq_NASDAQ - nu_NASDAQ^q * Kq)
est_N225 = cbind(q,bq_N225, nu_N225 ^q * Kq, bq_N225 - nu_N225^q * Kq)

b_estimates = round(rbind(est_SPX, est_DAX, est_NASDAQ, est_N225), digits=3)
write.table(b_estimates,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\b estimates.dat",
            sep=" & ", row.name=FALSE, col.names=FALSE)
