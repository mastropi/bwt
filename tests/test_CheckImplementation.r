# Created:      04-Nov-2015
# Author:       Daniel Mastropietro
# Description:  Check the BWT program as the implementation progresses. The first instance of these tests
#               is the change into a matrix-based computations.
#

library(waveslim)
library(wavethresh)

################### TEST 1: 04-Nov-2015 - CHECK MATRIX-BASED IMPLEMENTATION -----------------
# Load results from latest valid run in 2011...
load("E:/Daniel/Projects/P01-BWT/BWTRProject/workspace/20111105_FixKDer&KDer2CalculationsV12b.RData")

### 1.- Run the V12b version ---------------------------------------------------------------------
# Note however that this version is represented by the bwtV12c.r program. In fact this program was created
# equal to bwtV12b.r with the addition of global variables created for comparison purposes with the
# new matrix-based implementation (these global variables are for instance: gK, gK_der, gK_der2)
source("E:/Daniel/Projects/P01-BWT/BWTRProject/bwtV12c.r")
x_bwtV12b = bwt(x, prior=prior, wf="la16", adaptive=adaptive, fixed=fixed, fixed.values=c(0.5,1), max.length=100, margin=5, bands=TRUE)
## 2015/11/04: OK, it runs without errors

# Plot the signal, estimation and probability bands
with(x_bwtV12b, {
  x_est <<- x_hat_med;
  x_low <<- x_hat_low[, "005"]
  x_upp <<- x_hat_upp[, "995"]
}
)
ylim = range(c(x, x_low, x_est, x_upp))
plot(x, type="l", ylim=ylim)
lines(x_est, col="red")           # Estimated signal
lines(x_upp, col="red", lty=2)    # Upper bound
lines(x_low, col="red", lty=2)    # Lower bound
lines(x0, col="blue")             # True signal

### 2.- Run the matrix-based version ----------------------------------------
# Currently (2015/11/04) this is implemented in test_ImplementMatrixOperationsAllCases.r
# which is not prepared to be run "blindly". What I do now is open that code and run the part that
# is solely part of the process, this means:
# Start at line: "# Make a copy of the variables..."
# End at line 244: "par(mfrow=c(8,8), mar=c(0,0,0.5,0))" EXCLUDED


### 3.- Compare results (1) with (2) ----------------------------------------
summary(as.vector(K - gK))
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-0.0007250 -0.0007250 -0.0007250 -0.0006344 -0.0005608 -0.0003741
summary(as.vector(K_der - gK_der))
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-1.189e-05 -7.480e-10  0.000e+00  1.693e-09  6.910e-10  9.921e-06
summary(as.vector(K_der2 - gK_der2))
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-7.146e-05 -2.920e-09  0.000e+00  8.520e-09  3.470e-09  4.396e-05

## 2015/11/04: OK, the differences between results (1) and (2) are very small!
################### TEST 1: 04-Nov-2015 - CHECK MATRIX-BASED IMPLEMENTATION -----------------



################### TEST 2: 26-Nov-2015 - CHECK MATRIX-BASED IMPLEMENTATION -----------------
# This test is done using all the functions defined in the package
# Compile the needed R code
files = list.files(path="./R", pattern=".r$")
for (f in files) { cat("Compiling ./R/", f, "\n", sep=""); source(file.path("./R",f)) }

x_bwt = bwt(x, prior=prior, wf="la16", adaptive=adaptive, fixed=fixed, fixed.values=c(0.5,1), max.length=100, margin=5, bands=TRUE)

with(x_bwt, {
  x_est <<- x_hat_med;
  x_low <<- x_hat_low[, "005"]
  x_upp <<- x_hat_upp[, "995"]
}
)
ylim = range(c(x, x_low, x_est, x_upp))
plot(x, type="l", ylim=ylim)
lines(x_est, col="green", lwd=2)          # Estimated signal
lines(x_upp, col="green", lty=2, lwd=2)   # Upper bound
lines(x_low, col="green", lty=2, lwd=2)   # Lower bound
lines(x0, col="blue")                     # True signal
################### TEST 2: 26-Nov-2015 - CHECK MATRIX-BASED IMPLEMENTATION -----------------


# Compare with wave_band() from the wavetresh package
library(waveband)
x_wb = wave.band(x)   # Note: the intervals dranw are 99% probability bands (see ?wave.band)
t = (1:length(x))/length(x)
lines(t, x0[1:length(x)], col="black")

# Add the estimates from my code
lines(t, x_est, col="cyan")                   # Estimated signal
lines(t, x_upp, col="cyan", lty=2, lwd=2)     # Upper bound
lines(t, x_low, col="cyan", lty=2, lwd=2)     # Lower bound


x_bwt_normal = bwt(x, prior="normal", wf="la16", adaptive=adaptive, fixed=fixed, fixed.values=c(0.5,1), max.length=100, margin=5, bands=TRUE)
x_bwt_laplace = bwt(x, prior="laplacian", wf="la16", adaptive=adaptive, fixed=fixed, fixed.values=c(0.5,1), max.length=100, margin=5, bands=TRUE)


color = "orange"
with(x_bwt_normal, {
  x_est <<- x_hat_med;
  x_low <<- x_hat_low[, "005"]
  x_upp <<- x_hat_upp[, "995"]
}
)
ylim = range(c(x, x_low, x_est, x_upp))
plot(x, type="l", ylim=ylim)
lines(x_est, col=color, lwd=2)          # Estimated signal
lines(x_upp, col=color, lty=2, lwd=2)   # Upper bound
lines(x_low, col=color, lty=2, lwd=2)   # Lower bound
lines(x0, col="blue")                   # True signal
