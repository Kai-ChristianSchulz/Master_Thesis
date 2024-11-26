# fBM_simulation.R

################################################################################
# Pre-computing and saving fBM from t=0 to t=5052 with d= 10 and 100 simulations
# per day and H = 0.129, corresponding to the estimated H value for SPX.
################################################################################
# Note that instead of FBM() from the fractionalBM package, the circulant 
# embedding method which is implemented in "fBM_simulation.R" can be adopted. 
# It is much faster than FBM(). Since we implemented it late into writing the 
# thesis, we still used FBM() at this stage.
################################################################################

# library(devtools)
# install_github("732jhy/fractionalBM")
library(fractionalBM) # Package for computing fBM

n = 5052
H = 0.129

# d = 1/10 #####################################################################
d = 1/10
N = 1/d * n

set.seed(123)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation2.txt")

set.seed(1308)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation3.txt")

# d = 1/100 ####################################################################
d = 1/100
N = 1/d * n

set.seed(123)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation4.txt")

set.seed(1308)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation5.txt")

set.seed(1999)
B = c(0,FBM(N,H,t=n))
write(B, "D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBB-Simulation.txt")

