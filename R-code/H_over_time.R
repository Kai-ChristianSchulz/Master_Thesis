# H_over_time.R

################################################################################
# We investigate how H changes over time. First by using the same linear regression 
# method as in "smoothness.R", but now with a half year (approx 126 trading days)
# rolling window instead. Second, by splitting the data set into three parts and
# then applying the linear regression method to each third.
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


################################################################################
# First, for the SPX index, we use the same linear regression 
# method as in "smoothness.R", but now with a half year (approx 126 trading days)
# rolling window instead. We also restrict ourself to lags Delta = 1,...,10.
# All in all this calculation method quite unstable and not recommended. In the
# thesis we refer to an alternative method used by Gatheral and Bennedsen.
################################################################################

# Initializing 
SPX = OM[OM$Symbol==".SPX",]
RV_SPX = sqrt( SPX$rv5[! SPX$rv5 %in% 0] )
q = c(0.5,1,1.5,2,3)
delta = seq(1:10)
H = vector()

# We choose a rolling window of size 126 days (half a year of trading days)
window_start = 1
window_end = 126

# Calculating H value via linear regression for the rolling window
while(window_end <= length(RV_SPX)){
  data = RV_SPX[window_start:window_end]
  output = matrix(nrow = length(q), ncol = length(delta))
  for (i in 1:length(q)) {
    for (j in 1:length(delta))
      output[i,j] =  m(q[i],delta[j],data)
  }
  aq = vector()
  for (i in 1:length(q))
    aq = append(aq,lm(log(output[i,])  ~ log(delta))$coefficients[2] )
  H = append(H,lm(aq ~ q)$coefficients[2])
  
  window_start = window_start + 1
  window_end = window_end +1
}


# Plotting
plot(1:length(H), H, type="l",xaxt = 'n', ylab="H", xlab= "", main="SPX")
axis(side=1, at = seq(1,length(H), by = 252) , labels =seq(2001,2020, by=1))

################################################################################
# Second, we split the time horizon into the three parts 2000-2006, 2007-2013 
# and 2014-2020 and estimate H for each third based on the linear regression 
# method used in "smoothness.R". We repeat this for all indices in the OM 
# dataset.
# Note that not all indices have the same number of RV5 entries, they can 
# slightly vary. However, individually determining the correct breaks for each 
# index would be too much effort. We proxy by just dividing the length of each
# index by 3.
################################################################################

# Initializing #################################################################
indices = unique(OM$Symbol)

# Matrix that stores the H values for each index rowwise and each third columnwise
H_all = matrix(ncol=3,nrow=length(indices)) 
rownames(H_all) = indices

q = c(0.5,1,1.5,2,3)
delta = seq(1:50)

# Calculating H estimates for each index #######################################
for(index in indices) {
  
  breaks = c(0,1,2,3) * floor(length(OM[OM$Symbol==index,]$rv5)/3)
  
  for ( k in 1:3){
    
  data = sqrt ( OM[OM$Symbol==index,]$rv5[(breaks[k]+1):breaks[k+1]] )
  data = data[! data %in% 0] 
  output = matrix(nrow = length(q), ncol = length(delta))
  for (i in 1:length(q)) {
    for (j in 1:length(delta))
      output[i,j] =  m(q[i],delta[j],data)
  }
  xi = vector()
  for (i in 1:length(q))
    xi = append(xi,lm(log(output[i,])  ~ log(delta))$coefficients[2] )
  
  H_all[index,k] = lm(xi ~ q)$coefficients[2]
  }
  
}

# Saving as a table
H_all
write.table(round(H_all,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\H over time.dat", 
            sep=" & ")

