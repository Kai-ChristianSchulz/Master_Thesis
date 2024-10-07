# Plotting fBM for H = 0.2 and H = 0.8

library(fractionalBM)

set.seed(1999)
B_0.2 = c(0,FBM(1000,0.2,1000))
B_0.8 = c(0,FBM(1000,0.8,1000))

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB_0.2.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.8)
plot(B_0.2, type="l", xlab="t", ylab="", main = "H = 0.2")
text(-190,mean(B_0.2), expression(B[t]^H), xpd = TRUE, cex = 1.5) # Manually add y-axis labeling
dev.off()

postscript(file = "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB_0.8.eps", width = 6, height = 6, paper = "special", horizontal = FALSE)
par(cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.8)
plot(B_0.8, type="l", xlab="t", ylab="", main = "H = 0.8")
text(-190,mean(B_0.8), expression(B[t]^H), xpd = TRUE, cex = 1.5) # Manually add y-axis labeling
dev.off()

