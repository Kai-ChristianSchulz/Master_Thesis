# smoothness_all_indices.R

################################################################################
# We estimate the smoothness parameter H as the slope of the linear regression
# a_q vs q from "smoothness.R". We do this for all indices in the OM data set
################################################################################

# Initializing

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

indices = unique(OM$Symbol)

H = double(length(indices))
intersec = double(length(indices))
names(H) = indices
names(intersec) = indices

q = c(0.5,1,1.5,2,3)
delta = seq(1:50)

# Calculating H for all indices
for(index in indices) {
  data = sqrt ( OM[OM$Symbol==index,]$rv5 )
  data = data[! data %in% 0] # Removing 0 entries, because log(0) = infty
  output = matrix(nrow = length(q), ncol = length(delta))
  for (i in 1:length(q)) {
    for (j in 1:length(delta))
      output[i,j] =  m(q[i],delta[j],data)
  }
  xi = vector()
  for (i in 1:length(q))
    xi = append(xi,lm(log(output[i,])  ~ log(delta))$coefficients[2] )
  H[index] = lm(xi ~ q)$coefficients[2]
  intersec[index] = lm(xi ~ q)$coefficients[1]
}

# Saving as table
H_estimates = round(cbind(intersec, H), digits=3)
write.table(H_estimates,"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\H estimates.dat", 
            sep=" & ")
