# H over time

# We estimate H over time, using a the same linear regressiom method as in 
# "Estimating Smoothness", but now with a half year (approx 126 trading days) 
# rolling window  instead. 

# We examine the SPX

# We restrict ourself to Delta = 1,...,10

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
q = c(0.5,1,1.5,2,3)
delta = seq(1:10)
H = vector()

# We choose a rolling window of size 126 days (half a year of trading days)
window_start = 1
window_end = 126

# Calculating H value via linear regression for each rolling window
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


# Splitting the data set into three parts and calculating H for each  
# index in each third

# For example for SPX, we would do
# Entries 1:1749 -> 2000-2006
# Entries 1750:3511 -> 2007-2013
# 3512:length(dates) -> 2014-Beginning 2020
# breaks = c(0,1749,3511,5052)
# However, not all Indices have the same number of RV5 entries, thus
# we need to individually calculate the breaks. That would be too much effort,
# we proxy by just dividing the length of each index by 3.

indices = unique(OM$Symbol)
# Matrix that stores the H values for each index rowwise and each third columnwise
H_all = matrix(ncol=3,nrow=length(indices)) 
rownames(H_all) = indices

q = c(0.5,1,1.5,2,3)
delta = seq(1:50)


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

H_all
write.table(round(H_all,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\H over time.dat", 
            sep=" & ")

