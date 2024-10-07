# library(devtools)
# install_github("732jhy/fractionalBM")

# Simulating a fBM from t=0 to t=5052 with 100 Simulation per Day
# All in all simulating those 500 000 points takes about 15 Minutes
library(fractionalBM)

set.seed(1999)

d = 1/100
n = 5052
N = 1/d * n
H = 0.129


B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation.txt")

# d = 1/10
set.seed(123)
d = 1/10
N = 1/d * n
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation2.txt")

set.seed(1308)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation3.txt")

# d = 1/100
set.seed(123)
d = 1/100
N = 1/d * n
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation4.txt")

set.seed(1308)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation5.txt")

