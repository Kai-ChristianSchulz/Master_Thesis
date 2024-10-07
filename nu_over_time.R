# Nu over time

# Splitting the data set into three parts and calculating Nu for each  
# index in each third

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = unique(OM$Symbol)

Nu_all = matrix(ncol=3,nrow=length(indices)) 
rownames(Nu_all) = indices

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

Nu_all
write.table(round(Nu_all,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Nu over time.dat", 
            sep=" & ")