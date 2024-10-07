# Estimating Nu for all indices

# For each index we estimate Nu as the empirical standard deveation of the 
# increments of the log-volatility process with lag Delta = 1

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")

indices = unique(OM$Symbol)
Nu = double(length(indices))
names(Nu) = indices

for (ind in indices){
  data = sqrt ( OM[OM$Symbol==ind,]$rv5 )
  data = data[! data %in% 0] 
  sigma_2 = data[seq(2,length(data), by=1)]
  sigma_1 = data[seq(1,length(data)-1, by=1)]
  increments = log(sigma_2) - log(sigma_1)
  Nu[ind] = sd(increments)
}

Nu

write.table(round(Nu,digits=3),"D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\Nu estimates.dat", 
            sep=" & ")
