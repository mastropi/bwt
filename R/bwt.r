#**********************************************************************************************
#*************************************** MAIN PROGRAM *****************************************

# DEPENDENCIES
# ============
# waveslim
# wavethresh
#

# INDEX
# =====
# BAYES WAVELET THRESHOLDING
# bwt


# TODO: (2011/08/31) Before starting doing any calculations on y or x, standardize these functions
# by dividing by sigma. In fact, all expressions are actually functions of y/sigma, and not of y
# and sigma separately! (e.g. B(u), etc.)
#


#################################### BAYES WAVELET THRESHOLDING ###############################
bwt = function(x, params=NULL, prior="normal", adaptive=FALSE, wf="la16", fixed=TRUE, fixed.values=c(0.5, 1), max.length=512, margin=20, nro.points=20, bands=TRUE, logfile=NULL, header=NULL)
## x:               The signal to estimate via bayesian wavelet thresholding.
## params:          Vector with the initial estimates for the parameters of the posterior distribution of theta,
##                  when adaptive = FALSE. Otherwise, this should be left NULL.
##                  The parameters are, in this order:
##                  - p:           Probability of delta(theta) function in normal mixture model
##                  - sigma:       Value of standard deviation of noise
##                  - tau:         Value of standard deviation of the prior for theta
##                  The default values of the parameters correspond to the normal prior.
## prior:           Prior distribution to use for theta ("normal", "laplacian")
## adaptive:        Whether the 2 hyperparameters (p, tau) should be estimated adaptively, one
##                  for each decomposition level according to the form 2^(-alpha*j), 2^(-beta*j).
## wf:              Wavelet Filter to use. Possible values are any valid values for parameter wf in dwt() function.
##                  The default ("la16") is Daubechies Least Asymmetric wavelet of length 16 (N = 8, L = 16)
##                  according to notation in her papers of 1988 and 1992.
##                  This supersedes the options filter.number and family that were used in previous versions of bwt.
## fixed:           Logical. Whether the parameters alpha and beta in the ADAPTIVE framework are held
##                  fixed to the values given in fixed.values respectively (default to alpha = 0.5 and beta = 1).
## fixed.values:    2-D vector with the fixed values to be used for the parameters alpha and beta in the ADAPTIVE framework.
## max.length:      Maximum length of signal that allows using the theoretical expression for K_der and K_der2.
##                  If the length of the signal is larger than max.length, then a spline is fit to K(u)
##					in order to compute K_der and K_der2.
## margin:          Number of u-steps to leave at each side of the u sequence used to compute the cumulant generating function
##                  and its derivatives when they are computed using a smoothing spline fit to K(u)
##                  This is done to avoid inconsistent values of K_der and K_der2 that result at the borders in the
##                  spline prediction process.
##                  If there are inconsistent values of the above, try increasing the value of margin.
##                  This value cannot be increased too much because overflow errors may occur (due to large values of u).
## nro.points:      Nro. of points of u where K(u) and its derivatives are evaluated for the computation of the
##                  probability bands. Default: 20
## bands:           Logical. Whether the probability bands are requested or not. Default: TRUE
## logfile:         File where the log of the process is stored (in double quotes).
## header:          Header to include in the log file.
{
	### Display the function call:
	cat("The function bwt() was called with the following parameters:\n")
	cat("params=", params, "\n")
	cat("prior=", prior, "\n")
	cat("adaptive=", adaptive, "\n")
	cat("wf=", wf, "\n")
	cat("fixed=", fixed, "\n")
	cat("fixed.values=", fixed.values, "\n")
	cat("max.length=", max.length, "\n")
	cat("margin=", margin, "\n")
	cat("nro.points=", nro.points, "\n")
	cat("bands=", bands, "\n")
	cat("logfile=", logfile, "\n")
	cat("header=", header, "\n")
	
	### Parsing input parameters
	if (prior == "normal")
	{
		# The spline.df.factor is set to 100% in the normal case because there is no need to smooth the fit
		# of (r, K_der) since the values of K_der and K_der2 are correctly estimated in the case of the normal
		# prior.
		spline.df.factor = 1;
	}
	
	# Log file
	if (!is.null(logfile))
		sink(logfile)
	cat(header)
	
	####################### Centering and scaling the signal ##################################
	# The signal is centered to have zero mean and scaled to unit variance so that the scaling coefficient is zero
	# and the range of variation of u is defined.
	x_mean = mean(x);
	x_sd = sd(x);
	x0 = x;
	x = (x0 - x_mean) / x_sd
	####################### Centering and scaling the signal ##################################
	
	################################### Wavelet decomposition #################################
	cat("Computing the wavelet decomposition of the signal...\n")
	nro_levels = floor(log2(length(x)))
	y = dwt(x, wf=wf, n.levels=nro_levels, boundary="periodic")
	## Note that when the boundary="reflection", the coefficient vectors at each level are double as long as the coefficient vectors when boundary="periodic".
	y$N = length(x);                        # Add the information on the length of the signal to the dwt object.
	yD = as.vector(unlist(y)[1:(y$N-1)])    # The first N-1 coefficients of unlist(y) are the wavelet coefficients,
	# ordered from finest to coarsest scale.
	## NOTE that yD must be calculated BEFORE the addition of a character variable in the list. Otherwise, yD
	## becomes character!!
	y$nlevels = nro_levels;                 # Add the information on the number of levels of decomposition to the dwt object.
	y$wf = attr(y, "wavelet")               # Add the wavelet name used for easier extraction.
	################################### Wavelet decomposition #################################
	
	
	################################### Parameter estimation ##################################
	params = params_estimate(y, params=NULL, prior=prior, adaptive=adaptive, fixed=fixed, fixed.values=fixed.values)
	################################### Parameter estimation ##################################
	
	
	##################################### Signal estimation ###################################
	# Signal estimation via the inverse wavelet transform of the posterior median of theta
	cat("Estimating the signal via the inverse wavelet transform of the posterior median of theta...\n")
	theta_hat = posterior_median(y, params, prior=prior, adaptive=adaptive)
	
	# Remove extra elements in theta_hat, o.w. the inverse wavelet transform is incorrectly computed.
	theta_hat$N = NULL;
	theta_hat$nlevels = NULL;
	theta_hat$wf = NULL;
	x_hat = idwt(theta_hat);
	##################################### Signal estimation ###################################
	
	
	################################ Probability Bands for the signal #########################
	x_hat_low = NULL;
	x_hat_med = NULL;
	x_hat_upp = NULL;
	if (bands)
	{
		x_bands = prob_bands(y, params, prior=prior, max.length=max.length, margin=margin, nro.points=nro.points)
		x_hat_low = x_bands$x_hat_low;
		x_hat_med = x_bands$x_hat_med;
		x_hat_upp = x_bands$x_hat_upp;
		
		# De-centering and re-scaling
		x_hat_low = x_hat_low*x_sd + x_mean;
		x_hat_med = x_hat_med*x_sd + x_mean;
		x_hat_upp = x_hat_upp*x_sd + x_mean;
	}
	################################ Probability Bands for the signal #########################
	
	
	#################### De-centering and re-scaling of reconstructed signal ##################
	x = x0;
	x_hat = x_hat*x_sd + x_mean;
	#################### De-centering and re-scaling of reconstructed signal ##################
	
	if (!is.null(logfile))
		sink()
	
	list(x=x, x_hat=x_hat, x_hat_low=x_hat_low, x_hat_med=x_hat_med, x_hat_upp=x_hat_upp, theta_hat=theta_hat, y=y, extension="periodic", params=params, prior=prior, logfile=logfile)
}
#################################### BAYES WAVELET THRESHOLDING ###############################


#*************************************** MAIN PROGRAM *****************************************
#**********************************************************************************************
