# fBM_circulant_embedding.R

################################################################################
# We implement the circulant embedding method for simulating fBM proposed 
# by Helgason, Pipiras and Abry in "Fast and exact synthesis of stationary 
# multivariate Gaussian time series using circulant embedding"-
# This method uses Fast Fourier Transform to compute N sample points of fBM
# in O(log(N)*N). 
################################################################################

################################################################################
# First we simulate the increments of the fBM (Gaussian noise) X and then compute
# fBM as cumulative sum of the increments. We follow the steps from Helgason et
# al. Note that we simulate fBM at the equidistant time points t=1,...,N.
################################################################################

# Chosen parameters
N = 500000
H = 0.129

# 0) Autocovarinace function for Gaussian noise
gam_tilde = function(tau,H){
  1/2 * (abs(tau+1)^(2*H) + abs(tau-1)^(2*H) - 2* abs(tau)^(2*H))
}

# 1) Covariance embedding
r = sapply(c(0:(N-1),N,(-N+1):-1), function(tau) gam_tilde(tau,H))

# 2) Transfer to spectral domain
lam = Re(fft(r))

# 3) Cholesky factorization
# Can be skipped since we are in the univariate case and simply A_m = sqrt(2)

# 4) Noise Generation
W0 = rnorm(2*N,mean=0, sd=sqrt(1/2) )
W1 = rnorm(2*N,mean=0, sd=sqrt(1/2) )
W = W0 + 1i*W1
Z = sqrt(2) * W

# 5) Spectral Synthesis
X = Re(fft(sqrt(lam)*Z) * 1/sqrt(2*N))[1:N] 

# Calculating B as cumsum, saving and (optionally) plotting
B_circ = cumsum(X)
write(B_circ, paste0("D:\\Aachen\\Mathe\\12.Semester\\Masterarbeit\\R-Programme\\fBM_H-",H,"_N-",N,".txt") )
# plot(B_circ, type="l", col="blue", ylab="")


################################################################################


# Checking results by comparing them to another simulation method
library(fractionalBM)
B = FBM(N,H,t=N)
ymin = min(path_B,B)
ymax = max(path_B,B)
plot(path_B, type="l", col="blue", ylab="", ylim=c(ymin,ymax), xlab="")
lines(B, col="red")