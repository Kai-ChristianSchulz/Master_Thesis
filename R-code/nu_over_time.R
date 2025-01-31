# nu_over_time.R

################################################################################
# We investigate how nu changes over time. To this end we split the time horizon
# into three parts 2000-2006, 2007-2013 and 2014-2020 and estimate nu for each
# third as the standard deviation of the increments with lag Delta = 1 of the 
# log-realized volatility. We repeat this for each index in the OM data set.
# Note that not all indices have the same number of RV5 entries, they can 
# slightly vary. However, individually determining the correct breaks for each 
# index would be too much effort. We proxy by just dividing the length of each
# index by 3.
################################################################################


# Initializing #################################################################
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = unique(OM$Symbol)

# Matrix that stores the nu estimates for each index rowwise and each third columnwise
Nu_all = matrix(ncol=3,nrow=length(indices)) 
rownames(Nu_all) = indices

# Calculating nu estimates
for(index in indices) {
  
  breaks = c(0,1,2,3) * floor(length(OM[OM$Symbol==index,]$rv5)/3)
  
  for ( k in 1:3){
    
    data = sqrt ( OM[OM$Symbol==index,]$rv5[(breaks[k]+1):breaks[k+1]] )
    data = data[! data %in% 0] 
    sigma_2 = data[seq(2,length(data), by=1)]
    sigma_1 = data[seq(1,length(data)-1, by=1)]
    increments = log(sigma_2) - log(sigma_1)
    Nu_all[index,k] = sd(increments)
  }
  
}

# Saving as table
Nu_all
write.table(round(Nu_all,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Nu over time.dat", 
            sep=" & ")
