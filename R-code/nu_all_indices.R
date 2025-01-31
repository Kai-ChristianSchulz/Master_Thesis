# nu_all_indices.R

################################################################################
# We estimate nu as the standard deviation of the increments with 
# lag Delta = 1 of the 5 minute log-realized volatility. We repeat this for 
# all indices in the OM data set.
################################################################################


# Initializing #################################################################
OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = unique(OM$Symbol)
Nu = double(length(indices))
names(Nu) = indices

# Calculating nu estimate for each index #######################################
for (ind in indices){
  data = sqrt ( OM[OM$Symbol==ind,]$rv5 )
  data = data[! data %in% 0] 
  sigma_2 = data[seq(2,length(data), by=1)]
  sigma_1 = data[seq(1,length(data)-1, by=1)]
  increments = log(sigma_2) - log(sigma_1)
  Nu[ind] = sd(increments)
}

Nu

# Saving as table
write.table(round(Nu,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Nu estimates.dat", 
            sep=" & ")
