# Estimating log [ E( sigma_t *sigma_(t+Delta) ) ] as a function of Delta^(2H)
# for indices SPX, DAX, NASDAQ, Nikkei225

OM = read.csv("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Daten\\oxfordmanrealizedvolatilityindices.csv")
indices = c(".SPX",".GDAXI",".IXIC", ".N225")

H = c(0.129,0.098,0.126,0.119)
names(H) = indices

vol = vector(mode="list", length = length(indices))
names(vol) = indices

# Extracting volatilitiy data
for (index in indices){
  data = sqrt ( OM[OM$Symbol==index,]$rv5 )
  data = data[! data %in% 0]
  vol[[index]] = data
}

Delta = 1:50
AutoE = matrix(nrow=length(Delta), ncol=length(indices))
colnames(AutoE) = indices 


for (index in indices){
  n = length(vol[[index]])
  for (del in Delta){
    AutoE[del,index] = 1/(n-del)* sum(vol[[index]][1:(n-del)] * vol[[index]][(1+del):n])
  }
}


# Plotting 
mains = c("SPX, H = 0.129 ","GDAXI, H = 0.098","IXIC, H = 0.126", "N225, H = 0.119")
names(mains) = indices
file_names = c("SPX","GDAXI","IXIC","N225")
names(file_names) = indices

for (index in indices){
  postscript(file = paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\Bilder\\AutoExp_",file_names[index],".eps"), width = 6, height = 6, paper = "special", horizontal = FALSE)
  par(cex.lab = 1.3 , cex.axis = 1.2, cex.main = 1.5,mar = c(5, 5, 4, 2) + 0.1)
  x = Delta ^(2*H[index])
  plot(x,log(AutoE[,index]),type="p", pch=16, xlab=expression(Delta^{2 *H}),  
       ylab = expression(log (E(sigma[t]*sigma[t + Delta]))), main=mains[index] )
  abline(lm(log(AutoE[,index]) ~ x), col="red")
  dev.off()
}


