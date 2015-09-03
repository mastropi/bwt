#**********************************************************************************************
#****************************** FUNCTIONS FOR SIGNAL ESTIMATION *******************************

# INDEX
# =====
# FMARGINAL
# hyM
# fyM
#
# MARGINAL LIKELIHOOD
# loglike
# loglike_adaptive
#
# PARAMETER ESTIMATES
# map_params
# params_estimate
#
# KILL WAVELET COMPONENTS
# kill
#
# POSTERIOR MEDIAN
# posterior_median


# TODO: (2015/09/03) Put loglike() and loglike_adaptive() into a single function
# by adding a new parameter called adaptive.
#
# TODO: (2012/01/12) Function fmarginal(): Add new parameters to the function to pass the computed
# ratios needed for the calculation of the fmarginal, namely: sigma/tau and p/tau.
# This should accelerate the process because this function is called several times in the loop
# on j in function loglike_adaptive().
#


######################################## FMARGINAL ############################################
# Computes the pdf normalizing function defined as h(y) in the paper affecting the denominator of g(x/y)
hyM = function(y, params, prior="normal")
## y:           Number or vector at which h(y) is evaluated.
## params:      Matrix having as many rows as the length of y and 3 columns,
##				each equal, respectively, to the hyperparameter vectors p, sigma and tau
##				appearing in the posterior distribution of theta
##              (theta = mean of distribution of wavelet coefficients y).
##              The parameters are, in this order:
##              - p:           Probability of deltaDirac(theta) function in normal mixture model for corresponding scale j
##              - sigma:       Value of standard deviation of noise for corresponding scale j
##              - tau:         Value of standard deviation of the prior for theta for the corresponding scale j
##				In general sigma will be a constant for all scales.
## prior:		Name of the distribution function used as prior for theta.
##              Possible values: "normal", "laplacian"
{
	hy = NULL
	
	if (prior == "laplacian")
	{
		# Parse parameters of posterior distribution used in the calculation of h(y)
		sigma = params[,2];
		tau = params[,3];
		ysigma = y/sigma;
		sigmatau = sigma/tau;
		ytau = y/tau;
		
		hy = exp(ytau + log(pnorm(-ysigma-sigmatau))) + exp(-ytau + log(pnorm(ysigma-sigmatau)));
	}
	
	return(hy)
}

# Computes marginal pdf of y, the wavelet coefficients
fyM = function(y, params, prior="normal", hy=NULL)
## y:           See definition given in function hyM().
## params:      See definition given in function hyM().
## prior:       See definition given in function hyM().
## hy:			When parameter prior = "laplacian", this is a number or vector the
##				same dimension as y containing the pdf normalizing function defined as h(y)
##				in the paper as the denominator of g(x/y)
##				Otherwise use NULL.
##				Default: NULL
{
	# Parse parameters of posterior distribution
	p = params[,1];
	sigma = params[,2];
	tau = params[,3];
	sigmatau = sigma/tau;
	
	if (prior == "normal")
	{
		a = 1/(1 + sigmatau^2);
		fmarginal = (1-p)*dnorm(y, sd=sigma) + p*dnorm(y, sd=tau/sqrt(a));
	}
	else if (prior == "laplacian")
	{
		fmarginal = (1-p)*dnorm(y, sd=sigma) + 0.5*p/tau*exp(0.5*sigmatau^2)*hy;
		## Note that dnorm(y, sd=sigma) is equivalent to 1/sigma * dnorm(y/sigma) which is what is written in our paper.
	}
	
	return(fmarginal)
}
######################################## FMARGINAL ############################################


################################## MARGINAL LOGLIKELIHOOD #####################################
# Definition of the marginal loglikelihood to be maximized (used when calling R function optim() in function params_estimate())
loglike = function(theta, y, prior="normal")
{
	params = map_params(theta)
	hy = hy(y, params, prior=prior)
	-sum(log(fmarginal(y, params, prior, hy)))   # The minus sign is used because the function is minimized, not maximized (with optim()).
}

# Definition of the marginal loglikelihood to be maximized, in which the hyperparameters are different
# for each level of decomposition. This is used when parameter adaptive = TRUE.
loglike_adaptive <- function(theta, y, prior="normal", fixed.values=c(0.5, 1))
## y:               Wavelet decomposition object obtained with dwt()
## fixed.values:    2-D vector with the values of parameters alpha and beta when they are held fixed.
{
	sigma = exp(theta[1]);
	C1 = exp(theta[2]);
	C2 = exp(theta[3]);
	if (length(theta) == 5)     # This happens when the parameters alpha and beta are NOT held fixed.
	{
		alpha = exp(theta[4]);
		beta = exp(theta[5]);
	}
	else
	{
		# Fix the values of alpha and beta in order to have convergence (at least in Matlab it was like this)
		alpha = fixed.values[1];
		beta = fixed.values[2];
	}
	
	logL = 0;
	for (j in 0:(y$nlevels-1))
	{
		yj = y[[y$nlevels-j]]
		tauj = C1*2^(-alpha*j)
		pj = min(C2*2^(-beta*j), 1)
		hyj = hy(yj, c(pj, sigma, tauj), prior=prior)
		logL = logL + sum(log(fmarginal(yj, c(pj, sigma, tauj), prior, hyj)))
	}
	
	return(-logL)   # The minus sign is used because the function is minimized, not maximized (with optim()).
}
################################## MARGINAL LOGLIKELIHOOD #####################################


################################## PARAMETER ESTIMATES ########################################
# Transforms the parameters via functions that map the real line into the valid space for the
# parameter values
map_params = function(theta)
{
	p = 1/(1+exp(-theta[1]))
	sigma = exp(theta[2])
	tau = exp(theta[3])
	
	params = vector("numeric", length=3);
	params[1] = p;
	params[2] = sigma;
	params[3] = tau;
	
	return(params)
}

# Estimates the hyper parameters
params_estimate = function(y, params=NULL, prior="normal", adaptive=FALSE, fixed=TRUE, fixed.values=c(0.5, 1))
## y:               A dwt object containing the wavelet coefficients plus the following additional attributes:
##                  $N:         The length of the original signal.
##                  $nlevels:   The number of levels of the wavelet decomposition
##                  $wf:        The value of parameter wf used in the dwt() function that computed the wavelet coefficients.
## params:          List containing the initial parameter estimates for sigma, tau, p.
##                  The elements of the list must be sigma, tau and p, each of which is a real number.
## prior:           Prior distribution to use for theta, the mean of the wavelet coefficients y. Possible values: "normal", "laplacian"
## adaptive:        Whether the 2 hyperparameters (tau, p) should be estimated adaptively, one
##                  for each decomposition level according to the form 2^(-alpha*j), 2^(-beta*j), respectively for each parameter.
## fixed:           Logical. Whether the parameters alpha and beta in the ADAPTIVE framework are held
##                  fixed to the values given in fixed.values respectively (default to alpha = 0.5 and beta = 1).
## fixed.values:    2-D vector with the fixed values to be used for the parameters alpha and beta in the ADAPTIVE framework.
{
	cat("Fitting parameters of posterior distribution of theta (p, sigma, tau)...\n")
	
	# Create the vector containing the wavelet coefficients
	# First eliminate the character elements present in y (o.w. yD takes character values)
	y$wf = NULL;
	yD = as.vector(unlist(y)[1:(y$N-1)])
	
	if (!is.null(params))
	{
		sigma0 = params$sigma;
		tau0 = params$tau;
		p0 = params$p;
	}
	else
	{
		if (prior == "normal")
		{
			p0 = 0.5;
			sigma0 = 2;
			tau0 = 50;
		}
		else if (prior == "laplacian")
		{
			p0 = 0.5;
			sigma0 = 1;
			tau0 = 1;
		}
	}
	if (!adaptive)
	{
		fit = optim(c(log(p0/(1-p0)), log(sigma0),log(tau0)), loglike, y=yD, prior=prior)
		params = map_params(fit$par)
		p = params[1];
		sigma = params[2];
		tau = params[3];
		# Create a vector of length nlevels so that the parameters can be accessed easily in both the adaptive and non-adaptive cases.
		p = rep(p, y$nlevels)
		tau = rep(tau, y$nlevels)
	}
	else
	{
		sigma0 = 1;
		C10 = 50;
		C20 = 1;
		if (fixed)
		{
			fit = optim(c(log(sigma0), log(C10), log(C20)), loglike_adaptive, y=y, prior=prior, fixed.values=fixed.values)
		}
		else
		{
			alpha0 = fixed.values[1];
			beta0 = fixed.values[2];
			fit = optim(c(log(sigma0), log(C10), log(C20), log(alpha0), log(beta0)), loglike_adaptive, y=y, prior=prior)
		}
		
		# Parameter estimates obtained
		params <- exp(fit$par)
		sigma = params[1]
		C1 = params[2]
		C2 = params[3]
		if (fixed)
		{
			params[4] = fixed.values[1];
			params[5] = fixed.values[2];
		}
		alpha = params[4];
		beta = params[5];
		# Create a vector of length nlevels so that the parameters can be accessed easily in both the adaptive and non-adaptive cases.
		tau = rep(NA, y$nlevels);
		p = rep(NA, y$nlevels);
		for (j in 0:(y$nlevels-1))
		{
			tau[j+1] = C1*2^(-alpha*j)
			p[j+1] = min(1, C2*2^(-beta*j))
		}
	}
	params = list(sigma=sigma, tau=tau, p=p);
	
	return(params)
}
################################## PARAMETER ESTIMATES ########################################


############################### KILL WAVELET COMPONENTS #######################################
# Kills some of the wavelet components in order to recover the original signal
kill = function(y, params, prior="normal")
## y:           A vector containing the wavelet coefficients to kill.
## params:      A vector with parameters p, sigma, tau estimated from the marginal pdf.
## prior:       Prior distribution used for tau (either "normal" or "laplacian")
##
## RETURNED VALUES
## theta_hat:   A vector containing the coefficients in vector y after killing.
{
	p = params[1];
	sigma = params[2];
	tau = params[3];
	ysigma = y/sigma;
	sigmatau = sigma/tau;
	sigma2tau = sigma^2/tau;
	
	# Evaluation of the marginal pdf
	hy = hy(y, params, prior);
	fy = fmarginal(y, params, prior, hy);
	
	if (prior == "normal")
	{
		a = 1/(1 + sigmatau^2)
		b = sigma*sqrt(a)
		py = p*dnorm(y, sd=tau/sqrt(a)) / fy;
		Gy0 = pnorm(-a*y/b);
	}
	else if (prior == "laplacian")
	{
		# OPT-2011-START: Removed calculation of hy as this is computed above by calling function hy()
#        hy = exp(y/tau)*pnorm(-(ysigma + sigmatau)) + exp(-y/tau)*(1 - pnorm(-(ysigma - sigmatau)));
#        hy = exp(y/tau)*pnorm(-ysigma-sigmatau) + exp(-y/tau)*pnorm(ysigma-sigmatau);
		# OPT-2011-END
#        py = bound( 0.5*p/tau * hy * exp(0.5*sigma2tau^2) / fy, 1, "min");
		py = 0.5*p/tau * hy * exp(0.5*sigma2tau^2) / fy;
		Gy0 = exp(y/tau)/hy * pnorm(-(y + sigma2tau)/sigma);
	}
	
	enter = vector("numeric", length=length(y));
	theta_hat = vector("numeric", length=length(y));
	for (i in 1:length(y))
	{
		if (py[i] < 0 | py[i] > 1)
			cat("ERROR: py out of range\n")
		else if (py[i] <= 0.5)
		{
			enter[i] = 1;
			theta_hat[i] = 0;
		}
		else if (py[i]*Gy0[i] > 0.5)
		{
			enter[i] = 2;
			if (prior == "normal")
			{
				theta_hat[i] = a*y[i] + b*qnorm(0.5/py[i], 0, 1);
			}
			else if (prior == "laplacian")
			{
				theta_hat[i] = y[i] + sigma2tau + sigma*qnorm(0.5*hy[i]*exp(-y[i]/tau)/py[i], 0, 1);
			}
		}
		else if (py[i]*Gy0[i] < 0.5)
		{
			if ((1 - py[i]) + py[i]*Gy0[i] >= 0.5)
			{
				enter[i] = 3;
				theta_hat[i] = 0;
			}
			else
			{
				enter[i] = 4;
				if (prior == "normal")
				{
					theta_hat[i] = a*y[i] + b*qnorm(1-0.5/py[i], 0, 1);
				}
				else if (prior == "laplacian")
				{
					theta_hat[i] = y[i] - sigma2tau + sigma*qnorm(1-0.5*hy[i]*exp(y[i]/tau)/py[i], 0, 1);
				}
			}
		}
	}
	
	return(theta_hat)
}
############################### KILL WAVELET COMPONENTS #######################################


##################################### POSTERIOR MEDIAN ########################################
# Computes the posterior median to estimate theta from the wavelet coefficients of the signal
posterior_median = function(y, params, prior="normal", adaptive=FALSE)
## y:			A dwt object containing the wavelet coefficients plus the following additional attributes:
##              $N:         The length of the original signal.
##              $nlevels:   The number of levels of the wavelet decomposition
##              $wf:        The value of parameter wf used in the dwt() function that computed the wavelet coefficients.
## params:      A list containing the estimated hyperparameters (sigma, tau, p) at each level of decomposition.
##              $sigma:     Single number with the estimated value of the noise standard deviation
##              $tau:       Vector of length y$nlevels, with the estimated value of the scale parameter of
##                          the wavelet coefficient's prior at each level of decomposition.
##              $p:         Vector of length y$nlevels, with the estimated value of the mixture probability
##                          at each level of decomposition.
## prior:       Prior considered for theta ("normal", "laplacian")
## adaptive:    Whether the hyperparameters (p, tau) were estimated adaptively, one
##              for each decomposition level according to the form 2^(-alpha*j), 2^(-beta*j).
##
## RETURNED VALUES
## theta_hat:   A dwt object containing the coefficients in y after killing.
{
	# Parse input parameters
	sigma = params$sigma;
	tau = params$tau;
	p = params$p;
	
	# Initialize theta_hat as a dwt object
	theta_hat = y;
	
	# Eliminate element wf which has a character value before unlisting y because o.w. the result
	# obtained by the unlisting is not numeric.
	y$wf = NULL;
	# Create the vector containing the wavelet coefficients
	yD = as.vector(unlist(y)[1:(y$N-1)])    # Note that vector unlist(y) is ordered from finest to coarsest scale.
	if (!adaptive)
	{
		cat("sigma = ", sigma, "tau = ", tau[1], "p = ", p[1], "\n")
		# Computation of theta_hat
#        ykilled = kill(yD, c(p[1], sigma, tau[1]), prior);     # For modwt() case.
		ykilled = kill(yD, c(p[1], sigma, tau[1]), prior);
		# Store the result of the coefficients after killing into the dwt object theta_hat
		for (j in 0:(y$nlevels-1))
		{
#            theta_hat[[j+1]] = ykilled[(y$N*j+1):(y$N*(j+1))]  # For modwt() case.
			theta_hat[[j+1]] = ykilled[(2^y$nlevels*(1 - 2^(-j)) + 1) : (2^y$nlevels*(1 - 2^(-(j+1))))] # Fill theta_hat from finest to coarsest scale.
		}
		# Set to 0 the scaling coefficients.
#        theta_hat[[y$nlevels+1]] = rep(0, y$N);                # For modwt() case.
		theta_hat[[y$nlevels+1]] = 0;
	}
	else        # The hyperparameters were estimated adaptively at each different scale of the wavelet decomposition.
	{
		cat("sigma = ", sigma, "\n")
		for (j in 0:(y$nlevels-1))
		{
			yj = y[[y$nlevels-j]]
			cat("j =", j, "tauj = ", tau[j+1], "pj = ", p[j+1], "\n")
			# Computation of theta_hat at level j
			theta_hat_j = kill(yj, c(p[j+1], sigma, tau[j+1]), prior);
			theta_hat[[y$nlevels-j]] = theta_hat_j;
		}
	}
	
	return(theta_hat)
}
##################################### POSTERIOR MEDIAN ########################################


#****************************** FUNCTIONS FOR SIGNAL ESTIMATION *******************************
#**********************************************************************************************
