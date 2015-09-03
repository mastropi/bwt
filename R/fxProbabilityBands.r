#**********************************************************************************************
#****************************** FUNCTIONS FOR PROBABILITY BANDS *******************************

# INDEX
# =====
# AUXILIARY FUNCTIONS NEEDED FOR CGF
# B
# B_der
# B_der2
#
# CUMULANT GENERATING FUNCTION
# cgf
# cgf_der
# cgf_der2
# 
# WAVELET FUNCTIONS
# wavelet_functions
#
# RANGE OF U
# urange
#
# PROBABILITY BANDS
# prob_bands


# TODO: (2011/08/27) Replace the calls to function BOUND() which I think takes quite a lot of time
#	because it uses the function which() to find values that are larger or smaller than
# the given bound. Instead use functions pmin() and pmax() which implement this by a call to vapply()
#


############################## AUXILIARY FUNCTIONS NEEDED FOR CGF #############################
# B(u) function defined in paper
B = function(u, y, sigma, tau, y1u, y2u)
{
	# OPT-2011-START: Changed '1 - pnorm(-x)' with its equivalent pnorm(x) (to avoid one substraction), and created variable sigma2
	# DM-2011/08/31-START: Replaced exp(x)*pnorm(z) with exp(x + log(pnorm(z))) in order to avoid overflow of exp(x);
	sigma2 = sigma^2
#    B = exp((sigma^2*u + y)/B = exp((sigma2*u + y)/tau)*pnorm(y1u) + exp(-(sigma2*u + y)/tau)*pnorm(-y2u)tau)*pnorm(y1u) + exp(-(sigma^2*u + y)/tau)*(1 - pnorm(y2u))
#    B = exp((sigma2*u + y)/tau)*pnorm(y1u) + exp(-(sigma2*u + y)/tau)*pnorm(-y2u)
	arg = (sigma2*u + y)/tau
	B = exp(arg + log(pnorm(y1u))) + exp(-arg + log(pnorm(-y2u)));
	## Note that there is no problem when pnorm() = 0 for the log() function because log(0) = -Inf and therefore
	## the result of the exp() function is 0, as it is exp(-Inf) = 0, as it should be.
	# DM-2011/08/31-END:
	# OPT-2011-END
	return(B)
}

# B'(u) function defined in paper
B_der = function(u, y, sigma, tau, y1u, y2u)
{
	sigma2 = sigma^2
	sigmatau = sigma/tau
	# OPT-2011-START: Changed '1 - pnorm(-x)' with its equivalent pnorm(x) (to avoid one substraction), and created variable sigma2
	# DM-2011/11/05-START: Fixed the calculation of B_der in order to avoid OVERFLOW by using the log() function.
	# Note that before using the log() I need to check that its argument is >= 0... This is achieved by using the apply() functions.
	# I also replaced the expression sigma/tau with sigmatau and the argument of the exponential with variable arg
	# to avoid repeated calculations that make computations inefficient.
#    B_der = sigma * ( exp((sigma^2*u + y)/tau) * (sigma/tau*pnorm(y1u) - dnorm(y1u)) - exp(-(sigma^2*u + y)/tau) * (sigma/tau*(1 - pnorm(y2u)) - dnorm(y2u)) );
#    B_der = sigma * ( exp((sigma2*u + y)/tau) * (sigma/tau*pnorm(y1u) - dnorm(y1u)) - exp(-(sigma2*u + y)/tau) * (sigma/tau*pnorm(-y2u) - dnorm(y2u)) );
	factor1 = apply(cbind(0, sigmatau*pnorm(y1u)  - dnorm(y1u)), 1, max)	# Lower bound to 0 (to avoid e.g. -1E-300)
	factor2 = apply(cbind(0, sigmatau*pnorm(-y2u) - dnorm(y2u)), 1, max)	# Lower bound to 0 (to avoid e.g. -1E-300)
	arg = (sigma2*u + y)/tau
	B_der = sigma * ( exp( arg + log(factor1) ) - exp(-arg + log(factor2)) )
	# DM-2011/11/05-END
	# OPT-2011-END
	return(B_der)
}

# B"(u) function defined in paper
B_der2 = function(u, y, sigma, tau, y1u, y2u)
{
	sigma2 = sigma^2
	sigmatau = sigma/tau
	sigma2tau = sigma2/tau
	dn1 = dnorm(y1u)
	dn2 = dnorm(y2u)
	# OPT-2011-START: Changed '1 - pnorm(-x)' with its equivalent pnorm(x) (to avoid one substraction), and created variable sigma2
	# DM-2011/11/05-START: Did an equivalent change to the one done on cfg_der() function
#    B_der2 = sigma * ( exp((sigma^2*u + y)/tau) * (sigma^2/tau*(sigma/tau*pnorm(y1u) - 2*dnorm(y1u)) - sigma*dnorm(y1u)*y1u) + exp(-(sigma^2*u + y)/tau) * (sigma^2/tau*(sigma/tau*(1-pnorm(y2u)) - 2*dnorm(y2u)) + sigma*dnorm(y2u)*y2u) );
#    B_der2 = sigma * ( exp((sigma2*u + y)/tau) * (sigma2/tau*(sigma/tau*pnorm(y1u) - 2*dnorm(y1u)) - sigma*dnorm(y1u)*y1u) + exp(-(sigma2*u + y)/tau) * (sigma2/tau*(sigma/tau*pnorm(-y2u) - 2*dnorm(y2u)) + sigma*dnorm(y2u)*y2u) );
	factor1 = apply(cbind(0, sigma2tau*(sigmatau*pnorm(y1u) - 2*dn1) - sigma*dn1*y1u), 1, max)	# Lower bound to 0 (to avoid e.g. -1E-300)
	factor2 = apply(cbind(0, sigma2tau*(sigmatau*pnorm(-y2u) - 2*dn2) + sigma*dn2*y2u), 1, max) # Lower bound to 0 (to avoid e.g. -1E-300)
	arg = (sigma2*u + y)/tau
	B_der2 = sigma * ( exp(arg + log(factor1)) + exp(-arg + log(factor2)) );
	# DM-2011/11/05-END
	# OPT-2011-END
	return(B_der2)
}
############################## AUXILIARY FUNCTIONS NEEDED FOR CGF #############################


################################## CUMULANT GENERATING FUNCTION ###############################
### Cumulant generating function
cgf = function(u, y, params, prior="normal")
## Input parameters:
## u:           Vector or matrix giving the points where the cumulant is evaluated
## y:           Number representing the given value of the wavelet coefficient for which the CGF is evaluated
## params:      Vector with parameters appearing in the posterior distribution of theta
##              (mean of distribution of wavelet coefficients y).
##              These are, in this order:
##              - p:           Probability of delta(theta) function in normal mixture model
##              - sigma:       Value of standard deviation of noise
##              - tau:         Value of standard deviation of the prior for theta
## prior:       Name of the distribution function used as prior for theta.
##              Possible values: "normal", "laplacian"
{
	# Parse parameters of posterior distribution
	p = params[1];
	sigma = params[2];
	tau = params[3];
	ysigma = y/sigma;
	sigmatau = sigma/tau;
	
	# Quantities required for the computation of the cumulant generating function
	if (prior == "normal")
	{
		a = 1/(1 + sigmatau^2);
		b = sigma*sqrt(a);
		fy = fmarginal(y, params, prior, NULL);
#        py = bound( p*dnorm(y, sd=tau/sqrt(a)) / fy, 1, "min" );
		py = p*dnorm(y, sd=tau/sqrt(a)) / fy;
		## The function min is to avoid the value of py going above 1 (which happened for extreme values of y!!)
		## And having py > 1 generates a negative value of 1 - py and possibly a negative value of mgf below
		## which is not possible (since the MGF must be non-negative).
	}
	else if (prior == "laplacian")
	{
		# OPT-2011-START: Replaced the expression for hy with a call to function hy()
#        hy = exp(y/tau)*pnorm(-(ysigma+sigmatau)) + exp(-y/tau)*(1 - pnorm(-(ysigma-sigmatau)));
		# OPT-2011-END
		hy = hy(y, params, prior);
		fy = fmarginal(y, params, prior, hy);
		# DM-2011/08/31-START: Replaced the product of the exponential and other things with the exponential of a sum in order to avoid OVERFLOW!
		# Note that I also (re-)add a bound of py with 1 because I got cases where py was slightly larger than 1 (such that 1-py = -4.4E-16, i.e. slightly negative!!))
		# When not using this bound I got NaN values in the computation of cgf[indneg] below.
#        py = bound( 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy, 1, "min" );
		#py = 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy;
		py = min(exp( log(0.5*p/tau) + log(hy) - log(fy) + 0.5*sigmatau^2 ), 1);	# min() AVOIDS having py > 1!
		# DM-2011/08/31-END
		
		# DM-2011/08/27-START: Replaced the ratio py/hy by the expression above which is equal to py/hy
		# The reason is: be able to bound the ratio by MAXDOUBLE, which is not possible when doing py/hy
		# because NaN may be produced (i.e. not only Inf), and min(MAXDOUBLE,NaN) = NaN (and not MAXDOUBLE).
		# DM-2011/08/31-START: Replaced exp(x)*... by exp(x + log(...))
		#consty = py/hy;
		#consty = apply(cbind(MAXDOUBLE, 0.5*p/tau * exp(0.5*sigmatau^2)/fy), 1, min);
		consty = apply(cbind(MAXDOUBLE, exp(log(0.5*p/tau) - log(fy) + 0.5*sigmatau^2)), 1, min);
# DM-2013/06/05-TEMP: Global variable to store the value of consty and compare it with my calculation in the matrix implementation
		gconsty <<- consty
		# DM-2011/08/31-END
		# DM-2011/08/27-END
		y1u = -(ysigma + sigmatau*(1+u*tau));
		y2u = -(ysigma - sigmatau*(1-u*tau));
	}
	
	# Cumulant generating function
	# NOTE: (2011/08/28) In order to avoid overflow of the exponentials present in the cgf, K(u)
	# (such as exp[(sigma*u)^2] and exp(y*u)), I divide the vector u into 2 pieces:
	# - the piece where y*u is positive
	# - the piece where y*u is negative
	# In this way the expression I use to compute the CGF in each piece does not involve exponentials with positive arguments
	# which may lead to overflow (such as exp(y*u) for y*u > 0 and exp(-y*u) for y*u < 0).
	# In addition, since y and u appear (ONLY) multiplying in the exponential, I consider 2 cases in order to divide the vector u:
	# - y >= 0
	# - y < 0
	cgf = u*0;  # This is to initialize the object CGF with the same size and type of u (either a matrix or a vector)
	# In this way I don't need to worry about checking whether u is a vector or a matrix.
	if (y >= 0)
	{
		indpos = which(u>=0);
		indneg = which(u<0);
	}
	else
	{
		indpos = which(u<=0);
		indneg = which(u>0);
	}
	
	if (prior == "normal")
	{
		cgf[indpos] = bound( (a*y*u[indpos] + 0.5*(b*u[indpos])^2) + log( bound( (1-py)*exp(-(a*y*u[indpos] + 0.5*(b*u[indpos])^2)) + py, MINDOUBLE, "max" ) ), MAXDOUBLE, "min" );
		cgf[indneg] = bound( 0.5*(b*u[indneg])^2 + log( bound( (1-py)*exp(-0.5*(b*u[indneg])^2) + py*exp(a*y*u[indneg]), MINDOUBLE, "max" ) ), MAXDOUBLE, "min" );
	}
	else if (prior == "laplacian")
	# The following expressions were corrected on 08/03/07, by adding terms of the type (sigma*u)^2
	{
		cgf[indpos] = 0.5*(sigma*u[indpos])^2 + y*u[indpos] + log( (1-py)*exp(-0.5*(sigma*u[indpos])^2 - y*u[indpos]) + consty*B(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) );
		cgf[indneg] = 0.5*(sigma*u[indneg])^2 + log( (1-py)*exp(-0.5*(sigma*u[indneg])^2) + consty*exp(y*u[indneg])*B(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) );
	}
	
	return(cgf)
}

### First derivative of the cumulant generating function
cgf_der = function(u, y, params, prior="normal")
## Input parameters:
## u:           Vector or matrix giving the points where the cumulant is evaluated
## y:           Number representing the given value of the wavelet coefficient for which the CGF is evaluated
## params:      Vector with parameters appearing in the posterior distribution of theta
##              (mean of distribution of wavelet coefficients y).
##              These are, in this order:
##              - p:           Probability of delta(theta) function in normal mixture model
##              - sigma:       Value of standard deviation of noise
##              - tau:         Value of standard deviation of the prior for theta
## prior:       Name of the distribution function used as prior for theta.
##              Possible values: "normal", "laplacian"
{
	# Parse parameters of posterior distribution
	p = params[1];
	sigma = params[2];
	tau = params[3];
	ysigma = y/sigma;
	sigmatau = sigma/tau;
	
	# Quantities required for the computation of the cumulant generating function
	if (prior == "normal")
	{
		a = 1/(1 + sigmatau^2);
		b = sigma*sqrt(a);
		fy = fmarginal(y, params, prior, NULL);
#        py = bound( p*dnorm(y, sd=tau/sqrt(a)) / fy, 1, "min" );
		py = p*dnorm(y, sd=tau/sqrt(a)) / fy;
		## The function min is to avoid the value of py going above 1 (which happened for extreme values of y!!)
		## And having py > 1 generates a negative value of 1 - py and possibly a negative value of mgf below
		## which is not possible (since the MGF must be non-negative).
	}
	else if (prior == "laplacian")
	{
		# OPT-2011-START: Replaced the expression for hy with a call to function hy()
#        hy = exp(y/tau)*pnorm(-(ysigma+sigmatau)) + exp(-y/tau)*(1 - pnorm(-(ysigma-sigmatau)));
		hy = hy(y, params, prior);
		# OPT-2011-END
		fy = fmarginal(y, params, prior, hy);
#        py = bound( 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy, 1, "min" );
		# DM-2011/11/05-START: Did the same change I did in cgf() function
#        py = 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy;
		py = min(exp( log(0.5*p/tau) + log(hy) - log(fy) + 0.5*sigmatau^2 ), 1);	# min() AVOIDS having py > 1!
		# DM-2011/11/05-END
		
		# DM-2011/11/05-START: Did the same change I did in cgf() function
#        consty = py/hy;
		consty = apply(cbind(MAXDOUBLE, exp(log(0.5*p/tau) - log(fy) + 0.5*sigmatau^2)), 1, min);
		# DM-2011/11/05-END
		y1u = -(ysigma + sigmatau*(1+u*tau));
		y2u = -(ysigma - sigmatau*(1-u*tau));
	}
	
	# Derivative of the cumulant generating function
	cgf_der = u*0;  # This is to initialize the object CGF with the same size and type of u (either a matrix or a vector)
	# In this way I don't need to worry about checking whether u is a vector or a matrix.
	if (y >= 0)
	{
		indpos = which(u>=0);
		indneg = which(u<0);
	}
	else
	{
		indpos = which(u<=0);
		indneg = which(u>0);
	}
	if (prior == "normal")
	{
		cgf_der[indpos] = py*(a*y + b^2*u[indpos]) / ( (1-py)*exp(-(a*y*u[indpos] + 0.5*(b*u[indpos])^2)) + py );
		cgf_der[indneg] = py*exp(a*y*u[indneg])*(a*y + b^2*u[indneg]) / ( (1-py)*exp(-0.5*(b*u[indneg])^2) + py*exp(a*y*u[indneg]) );
	}
	else if (prior == "laplacian")
	# The following expressions were corrected on 08/03/07, by adding terms of the type (sigma*u)^2 in the cumulant
	{
		cgf_der[indpos] = consty * ( B(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos])*(y + sigma^2*u[indpos]) + B_der(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) ) / ( (1-py)*exp(-0.5*(sigma*u[indpos])^2 - y*u[indpos]) + consty*B(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) )
		cgf_der[indneg] = consty * exp(y*u[indneg]) * ( B(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg])*(y + sigma^2*u[indneg]) + B_der(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) ) / ( (1-py)*exp(-0.5*(sigma*u[indneg])^2) + consty*exp(y*u[indneg])*B(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) );
	}
	
	return(cgf_der)
}

### Second derivative of the cumulant generating function
cgf_der2 = function(u, y, params, K_der, prior="normal")
## Input parameters:
## u:           Vector or matrix giving the points where the cumulant is evaluated
## y:           Number representing the given value of the wavelet coefficient for which the CGF is evaluated
## params:      Vector with parameters appearing in the posterior distribution of theta
##              (mean of distribution of wavelet coefficients y).
##              These are, in this order:
##              - p:           Probability of delta(theta) function in normal mixture model
##              - sigma:       Value of standard deviation of noise
##              - tau:         Value of standard deviation of the prior for theta
## K_der:       Vector containing the value of the first derivative of the cumulant generating function
##              for each value in the vector u.
## prior:       Name of the distribution function used as prior for theta.
##              Possible values: "normal", "laplacian"
{
	# Parse parameters of posterior distribution
	p = params[1];
	sigma = params[2];
	tau = params[3];
	ysigma = y/sigma;
	sigmatau = sigma/tau;
	
	# Quantities required for the computation of the cumulant generating function
	if (prior == "normal")
	{
		a = 1/(1 + sigmatau^2);
		b = sigma*sqrt(a);
		fy = fmarginal(y, params, prior, NULL);
#        py = bound( p*dnorm(y, sd=tau/sqrt(a)) / fy, 1, "min" );
		py = p*dnorm(y, sd=tau/sqrt(a)) / fy;
		## The function min is to avoid the value of py going above 1 (which happened for extreme values of y!!)
		## And having py > 1 generates a negative value of 1 - py and possibly a negative value of mgf below
		## which is not possible (since the MGF must be non-negative).
	}
	else if (prior == "laplacian")
	{
		# OPT-2011-START: Replaced the expression for hy with a call to function hy()
#        hy = exp(y/tau)*pnorm(-(ysigma+sigmatau)) + exp(-y/tau)*(1 - pnorm(-(ysigma-sigmatau)));
		hy = hy(y, params, prior);
		# OPT-2011-END
		fy = fmarginal(y, params, prior, hy);
#        py = bound( 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy, 1, "min" );
		# DM-2011/11/05-START: Did the same change I did in cgf() function
#        py = 0.5*p/tau * hy * exp(0.5*sigmatau^2)/fy;
		py = min(exp( log(0.5*p/tau) + log(hy) - log(fy) + 0.5*sigmatau^2 ), 1);	# min() AVOIDS having py > 1!
		# DM-2011/11/05-END
		
		# DM-2011/11/05-START: Did the same change I did in cgf() function
#        consty = py/hy;
		consty = apply(cbind(MAXDOUBLE, exp(log(0.5*p/tau) - log(fy) + 0.5*sigmatau^2)), 1, min);
		# DM-2011/11/05-END
		y1u = -(ysigma + sigmatau*(1+u*tau));
		y2u = -(ysigma - sigmatau*(1-u*tau));
	}
	
	# Second derivative of cumulant generating function
	cgf_der2 = u*0; # This is to initialize the object CGF with the same size and type of u (either a matrix or a vector)
	# In this way I don't need to worry about checking whether u is a vector or a matrix.
	if (y >= 0)
	{
		indpos = which(u>=0);
		indneg = which(u<0);
	}
	else
	{
		indpos = which(u<=0);
		indneg = which(u>0);
	}
	if (prior == "normal")
	{
		cgf_der2[indpos] = py*((a*y + b^2*u[indpos])^2 + b^2) / ( (1-py)*exp(-(a*y*u[indpos] + 0.5*(b*u[indpos])^2)) + py ) - (K_der[indpos])^2;
		cgf_der2[indneg] = py*exp(a*y*u[indneg])*((a*y + b^2*u[indneg])^2 + b^2) / ( (1-py)*exp(-0.5*(b*u[indneg])^2) + py*exp(a*y*u[indneg]) ) - (K_der[indneg])^2;
	}
	else if (prior == "laplacian")
	{
		cgf_der2[indpos] = K_der[indpos]*(y + sigma^2*u[indpos] - K_der[indpos]) + consty * ( sigma^2*B(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) + (y + sigma^2*u[indpos])*B_der(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) + B_der2(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) ) / ( (1-py)*exp(-0.5*(sigma*u[indpos])^2 - y*u[indpos]) + consty*B(u[indpos], y, sigma, tau, y1u[indpos], y2u[indpos]) )
		cgf_der2[indneg] = K_der[indneg]*(y + sigma^2*u[indneg] - K_der[indneg]) + consty * exp(y*u[indneg]) * ( sigma^2*B(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) + (y + sigma^2*u[indneg])*B_der(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) + B_der2(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) ) / ( (1-py)*exp(-0.5*(sigma*u[indneg])^2) + consty*exp(y*u[indneg])*B(u[indneg], y, sigma, tau, y1u[indneg], y2u[indneg]) );
	}
	
	return(cgf_der2)
}
################################## CUMULANT GENERATING FUNCTION ###############################


###################################### WAVELET FUNCTIONS ######################################
### Computes wavelet functions
wavelet_functions = function(N, nlevels, wf)
## Input parameters:
## N:               Length of the signal being wavelet transformed.
## nlevels:         Number of levels used in the wavelet decomposition (value of parameter nlevels in dwt())
## wf:              Name of the wavelet used in the decomposition (value of parameter wf in dwt())
##
## RETURNED VALUES
## psi:             Matrix with y$N number of rows and y$N-1 columns.
##                  Each row is a time point where the original signal is defined and each column
##                  is indexed by a scale j and a shift k as follows: col = 2^j + k
##                  This means that each column col contains a wavelet function Psi(j,k), where
##                  j = floor(log2(col)) and k = col - j
##
## METHOD
## Each wavelet function Psi(j,k), where j is the scale and k is its location, is obtained by performing
## the IDWT to a synthesized wavelet transform whose d(j,k) = 1 and the rest are set to 0.
## This is done using the dwt() function (regardless of whether the translation invariant wavelet tranform
## is used to compute the input parameter y, because what we need to compute here is the wavelet functions
## and the described method for computing the wavelet functions is based on the traditional wavelet transform.
{
	# Initialize a dwt object with all coefficients equal to 0.
	theta = dwt(rep(0, N), wf=wf, n.levels=nlevels, boundary="periodic");
	
	# Put a 1 on each position of detail coefficients and reconstruct the signal corresponding to that
	# distribution of coefficients (which gives the desired wavelet function).
	psi = matrix(0, nrow=N, ncol=N-1);
	for (j in 0:(nlevels-1))
	{
		cat("j=", j, "\n")
		for (k in 0:(2^j-1))
		{
			theta[[nlevels-j]][k+1] = 1;
			psi[,2^j+k] = idwt(theta);
			## This way of indexing Psi indicates that the column 2^j+k of Psi corresponds
			## to the values of the wavelet psi(j,k) corresponding to level j and shift k.
			# Restore the zero in the modified coefficient
			theta[[nlevels-j]][k+1] = 0;
		}
	}
	
	return(psi)
}
###################################### WAVELET FUNCTIONS ######################################


######################################### RANGE OF U ##########################################
# Defines the values of u at which the cumulant K(u) and its derivatives are computed.
urange = function(umax, ustep, prior, margin=10)
{
	if (prior == "normal")
	{
		u = seq(-umax, umax, ustep);
		# Indices that indicate the first and last values of u considered in the computation of K(u).
		ufirst = 1;
		ulast = length(u);
		uind = ufirst:ulast;
	}
	else if (prior == "laplacian")
	{
		u = seq(-(umax+margin*ustep), umax+margin*ustep, ustep);
		## Note that the maximum values considered for u are augmented by margin*ustep. This is because the derivatives K_der and K_der2 are
		## not computed using the closed form, but they are estimated by prediction from a smoothing spline fit to K(u). It was seen that
		## the prediction of the derivatives is not good at the borders.
		# Indices that indicate the first value of u to be considered in the evaluation of K_der and K_der2 in order to avoid the border effects.
		ufirst = which(abs(u+umax)<1E-6);   # This is tantamount to finding the indices where u = -umax
		ulast = which(abs(u-umax)<1E-6);    # This is tantamount to finding the indices where u = umax
		uind = ufirst:ulast;
	}
	
	list(u=u, uind=uind);
}
###################################### RANGE OF U #############################################


#################################### PROBABILITY BANDS ########################################
# Computes the probability bands for the true signal
prob_bands = function(y, params, prior="normal", max.length=512, margin=20, nro.points=20)
## y:               A dwt object containing the wavelet coefficients plus the following additional attributes:
##                  $N:         The length of the original signal.
##                  $nlevels:   The number of levels of the wavelet decomposition
##                  $wf:        The value of parameter wf used in the dwt() function that computed the wavelet coefficients.
## params:          A list containing the estimated hyperparameters (sigma, tau, p) at each level of decomposition.
##                  $sigma:     Single number with the estimated value of the noise standard deviation
##                  $tau:       Vector of length y$nlevels, with the estimated value of the scale parameter of
##                              the wavelet coefficient's prior at each level of decomposition.
##                  $p:         Vector of length y$nlevels, with the estimated value of the mixture probability
##                              at each level of decomposition.
## prior:           Prior considered for theta ("normal", "laplacian")
## max.length:      See documentation for bwt()
## margin:          See documentation for bwt()
## nro.points:      See documentation for bwt()
##
## RETURNED VALUES
## x_hat_low
## x_hat_med
## x_hat_upp
{
	### Parse input parameters
	# Parameter estimates
	sigma = params$sigma;
	tau = params$tau;
	p = params$p;
	# Define whether a spline should be fit to K(u) in order to estimate K_der and K_der2.
	spline.fit = TRUE;
	if (prior == "normal" | y$N <= max.length)
	{
		spline.fit = FALSE;
	}
	
	# Wavelet functions
	cat("Computing the confidence bands of the signal...\n")
	if (spline.fit) {
		cat("Using spline fit for the derivatives of the cumulant generating function, K(u)\n")
	}
	else {
		cat("WITH NO spline interpolation for the derivatives of K(u)\n")
	}
	cat("\tWavelet functions...\n")
	psi = wavelet_functions(y$N, y$nlevels, y$wf);
	
	##### DELETE ######
#    # For each time ti, which (j,k) have Psi(j,k)(ti) equal to 0.
#    ind0 = which(psi==0, arr.ind=TRUE);
#    # Indices where the wavelet is not 0.
#    indnot0 = which(psi!=0, arr.ind=TRUE);
	
#    indnot0 = matrix(FALSE, nrow=y$N, ncol=ncol(psi))
#    for (i in 1:y$N)
#    {
#        indnot0[i,] = psi[i,] != 0;
#    }
	##### DELETE ######
	
	### Determine the range of u for which K(u) is computed.
	K_der2_0 = rep(0, y$N);
	for (j in 0:(y$nlevels-1))
	{
		for (k in 0:(2^j-1))
		{
			# Time indices where the wavelet function is non-zero.
			# Contributions to K"(u) come only from time points where the wavelet function psi(j,k,t) is not zero as the square
			# of psi is multiplying each term of the summation expression giving K"(u) (see paper p.485).
			indnot0 = which(psi[,2^j+k]!=0);	# Recall that each row is a time point, so this selects all the TIME points where psi(j,k,t) is NOT zero.
			# Value of y(j,k)
			yjk = y[[y$nlevels-j]][k+1];
			# The following are NUMBERS (not vectors).
			Kjk_der_0 = cgf_der(0, yjk, c(p[j+1], sigma, tau[j+1]), prior="normal");
			Kjk_der2_0 = cgf_der2(0, yjk, c(p[j+1], sigma, tau[j+1]), Kjk_der_0, prior="normal");
			K_der2_0[indnot0] = K_der2_0[indnot0] + psi[indnot0,2^j+k]^2 * Kjk_der2_0;
		}
	}
	
	# Range for u
	# These are the points where the cumulant generating function K(u) and its first and second derivatives are evaluated
	# for the computation of the probability bands.
	# Note that there is a difference between the normal prior case and the laplacian prior case, because in the latter, the
	# derivatives of K(u) are computed via an estimation based on the smoothing spline fit to K(u) and there are border effects
	# that are avoided as much as possible by extending the range of computation of K(u).
	# DM-2013/06/04: Added the MAX() function in the computation of umax in order to have one single u vector in order to carry out the matrix implementation
	# which requires a KM matrix where the rows represent the y axis and the columns represent the u axis, whose points have to be the same
	# for all y values (i.e. all i indices).
	umax = max(round(3.5/sqrt(K_der2_0)))    # The round() function is to have always u=0 as one of the evaluating points and thus avoid problems in Fu
	# when u is very close to 0 generating the possibility of v1 and v2 not being exactly 0 when they should be.
	if (nro.points < 10)
	{
		nro.points = 20;
	}
	ustep = 2*umax/ceiling(nro.points);
	
	# DM-2013/06/04-START: Moved the computation of vector u from the beginning of the FOR i loop below to here
	# as for the matrix implementation we need a single u vector for all y values (i.e. all i indices)
	# (in V12b version, the u vector was computed for each i using the umax value computed above for that i)
	# In addition I renamed u with uext in order to avoid overwriting u inside the i loop which messed everything up,
	# and because it is more convenient to call the actual u values as u and not as usubset (for instance).
	UVar = urange(umax, ustep, prior, margin);
	uext = UVar$u;
	uind = UVar$uind;
	u = uext[uind];
	print(uext)
	# DM-2013/06/04-END
	
	### Cumulant generating function
	# Initializing the matrices that will contain the median, and the lower and upper probability limits.
	x_hat_med = vector("numeric", length=y$N);
	x_hat_low = as.data.frame(matrix(nrow=y$N, ncol=3));
	x_hat_upp = as.data.frame(matrix(nrow=y$N, ncol=3));
	colnames(x_hat_low) = c("005", "025", "050");
	colnames(x_hat_upp) = c("950", "975", "995");
	
	# Construct yvec so that the wavelet coefficients are accessed directly through a linear index which grows as 2^j+k,
	# as is the case with columns in psi.
	yvec = NULL;
	for (j in 0:(y$nlevels-1))
	{
		yvec = c(yvec, y[[y$nlevels-j]])
	}
	
	# DM-2013/06/04-TEMP: Create matrix gK and its derivatives (g stands for 'global') in the global environment in order to store the values of K for all i's and u's
	# and compare with their values with the results from the matrix implementation started in Jan-2012 (see tasks.txt and test_ImplementationMatrixOperationsAllCases.r)
	# Note that the number of columns of gK* is given by the length of u (not uext which is the length of K above) 
	# which defines the actual range of u where the K(u) is evaluated (to avoid the border effects in the Laplacian prior case)
	gK <<- matrix(data=0, nrow=y$N, ncol=length(u))
	gK_der <<- matrix(data=0, nrow=y$N, ncol=length(u))
	gK_der2 <<- matrix(data=0, nrow=y$N, ncol=length(u))
	
	cat("\tCumulant generating function and its derivatives...\n")
# Counter of wrong cases with u*K'(u) - K(u) < 0 and of K"(u) < 0, which should NOT happen
	wrong_count = 0
	wrong_min_value = 0
	wrong_values = NULL
	K_der2_neg_count = 0
	K_der2_min_value = Inf
	K_der2_neg_values = NULL
	for (i in 1:y$N)		# i is the index for time
	{
		# Only go over the indices j and k that have contribution to K(u). These are given by those (j,k) for which Psi(j,k) is non-zero.
		indnot0 = which(psi[i,]!=0);	# This selects all the scales j and shifts k of psi(j,k,i) that are NOT zero for the current analyzed TIME i.
		# Above, where K"(0) is computed, indnot0 INSTEAD contained all the TIMES (instead of scales and shifts as is the case here)
		# where psi(j,k,t) are not zero for the corresponding analyzed scale j and shift k.
		cat("K: i =", i, "\n")
#print(indnot0)
#print(length(indnot0))
		# First and second derivative
		K = rep(0, length(uext));
		K_der = rep(0, length(uext));
		K_der2 = rep(0, length(uext));
		
		for (l in indnot0)		# l is the index that puts together both scale j and shift k into one signle index, representing the columns of the wavelet functions psi.
		{
			j = floor(log2(l));
			k = l - 2^j;
			yjk = yvec[l];
#            cat("\tl=", l, "j=", j, "k=", k, "\n")
			
			# Cumulant
			Kjk = cgf(psi[i,l]*uext, yjk, params=c(p[j+1], sigma, tau[j+1]), prior=prior);
			# Completing the matrix cgfjk with the cases when psi = 0.
			K = K + Kjk;
			
			# First and second derivative of cumulant
			if (!spline.fit)
			{
				Kjk_der = cgf_der(psi[i,l]*uext, yjk, c(p[j+1], sigma, tau[j+1]), prior);
				Kjk_der2 = cgf_der2(psi[i,l]*uext, yjk, c(p[j+1], sigma, tau[j+1]), Kjk_der, prior);
#DM-2011/11/05: Used to check the calculation of K_der and K_der2 done by --essentially-- the implementation in B_der() and B_der2() which are called by cgf_der() and cgf_der2()
#cat("\t\tKjk_der=", Kjk_der, "\n")
#cat("\t\tKjk_der2=", Kjk_der2, "\n")
				K_der = K_der + psi[i,l] * Kjk_der;
				K_der2 = K_der2 + (psi[i,l])^2 * Kjk_der2;
			}
		}
#            else
#           {
		# Here was the spline fit to K(u) and the computation of K_der and K_der2 via the derivative of the spline fit
		# which is now located below, after the iteration on l finishes...!
#            }
#        }
		if (spline.fit)
		{
# DM-2013/04/06: Check if there are any Inf or NaN in K(u)
			print(K)
			# Fit a spline to the cumulant generating function in order to compute its derivatives
			# Here we see the justification of using uext, so that the spline fit does not generate border effects
			# in the actual and final range of u used in computations.
			K_spline = smooth.spline(uext, K);
			# Estimate first derivative K_der
			K_der_pred = predict(K_spline, uext, deriv=1);
			K_der = K_der_pred$y;
			# Fit a spline to the first derivative
			K_der_spline = smooth.spline(uext, K_der);
			# Estimate seconde derivative K_der2
			K_der2_pred = predict(K_der_spline, uext, deriv=1);
			K_der2 = K_der2_pred$y;
			
			# Restrict the values of uext to the useful range given by u which avoids inconsistent values in K_der and K_der2.
			u = uext[uind]
			K = K[uind];
			K_der = K_der[uind];
			K_der2 = K_der2[uind];
		}
		
#print(cbind(u,K))
#if (!(j == 6 && k == 3 && l == 67 && i == 1)) {
		# Store variables in the Global Environment for debugging (could also use the <<- operator)
		psi <<- psi;
		assign("y", y, envir=.GlobalEnv);
		assign("yvec" , yvec , envir=.GlobalEnv);
		assign("sigma", sigma, envir=.GlobalEnv);
		assign("tau", tau, envir=.GlobalEnv);
		assign("p", p, envir=.GlobalEnv);
		assign("u" , u , envir=.GlobalEnv);
		uext <<- uext
		uind <<- uind
		assign("uK", cbind(u,K), envir=.GlobalEnv)
		assign("i", i, envir=.GlobalEnv)
		assign("j", j, envir=.GlobalEnv)
		assign("k", k, envir=.GlobalEnv)
		assign("l", l, envir=.GlobalEnv)
		
		# DM-2013/06/04-TEMP: Update the matrices K, K_der and K_der2
		gK[i,] <<- K
		gK_der[i,] <<- K_der
		gK_der2[i,] <<- K_der2
#}
		
#keypressed = readline("********** Press ENTER to continue **********")
		
# Check whether the saddlepoint approximation can be computed without error, and the number of cases where
# v1 = 0, which prevents the computation of the variable r which is ~ 1/v1.
		wrong = vector(length=nrow(psi), "numeric");
		zeros = vector(length=nrow(psi), "numeric");
		diff_min = Inf;     # Variable that stores the minimum difference between u*K_der and K, so that we know how
		# negative is that difference going to be in case there are negative values.
		diff = u*K_der - K
		indwrong = which(diff < 0);
		wrong = length(indwrong);
		diff_min = min(diff);
		zeros = length(which(u*K_der - K == 0));
		K_der2_min_value = min(K_der2_min_value, min(K_der2))
		if (sum(wrong) > 0) {
			wrong_count = wrong_count + 1
			wrong_min_value = min(wrong_min_value, diff_min)
			wrong_values = c(wrong_values, diff_min)
			cat("u*K_der - K < 0: ", sum(wrong), ", smaller value: ", diff_min, "\n")
			cat("u*K_der - K = 0: ", sum(zeros), "\n")
			cat("K_der2 < 0", length(which(K_der2 < 0)), ", smaller value: ", min(K_der2), "\n")
			if (length(which(K_der2<0)) > 0) {
				cat("\tWrong cases where K_der2 < 0\n")
				print(K_der2[K_der2 < 0])
				K_der2_neg_count = K_der2_neg_count + 1
				K_der2_neg_values = c(K_der2_neg_values, min(K_der2))
			}
		}
		# Saddlepoint approximation.
		cat("\tSaddlepoint approximation to the distribution of the signal given the wavelet coefficients...\n")
		v1 = sign(u)*sqrt( bound(2*(u*K_der - K), 0, "max") );
		v2 = u*sqrt( bound(K_der2, 0, "max") );
		r = v1 + 1/v1 * log(v2/v1);
		Fu = pnorm(r);     # There are as many NA's as the length of the signal which correspond to the cases u = 0.
		
		### Computing the lower and upper bounds of the CI for the signal
		### The method is based on the prediction of the quantiles -1.96, 0 and 1.96 given by a smoothing spline
		### fitted to the points (r(u), K'(u)).
		# Excluding NA's and Inf of r to avoid errors in the smoothing spline fit.
		ind_exclude = which(is.na(r) | abs(r)==Inf);
		if (length(ind_exclude) > 0)
		{
			rk.spline = try(smooth.spline(r[-ind_exclude], K_der[-ind_exclude]))#, df=spline.df.factor*length(r[i,-ind0])))
		}
		else
		{
			rk.spline = try(smooth.spline(r, K_der))#, df=spline.df.factor*length(r[i,])))
			## Note that if I use the -ind_exclude notation as above when ind_exclude is empty, the vector r[-ind_exclude] is also EMPTY! Why? I don't know.
		}
		if (inherits(rk.spline, "try-error"))
		{
			cat("ERROR: The smoothing spline could not be fit. The case will be skipped (i =", i, ")\n")
#print(ind_exclude)
#cat("Follow:\t\t v1, v2, r, K, K_der, K_der2, u*K_der - K\n")
#print(cbind(v1, v2, r, K, K_der, K_der2, u*K_der - K))
		}
		else
		{
			# (2011/11/05) The following numbers are the normal quantiles for 0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995
			# and they were hardcoded in order to reduce computation time...
#            quantiles = predict(rk.spline, c(qnorm(0.005), qnorm(0.025), qnorm(0.05), 0, qnorm(0.95), qnorm(0.975), qnorm(0.995)))  # The 0 in the middle is qnorm(0.5)
			q1 = -2.575829;		# quantile 0.005
			q2 = -1.959964;		# quantile 0.025
			q3 = -1.644854;		# quantile 0.05
			q4 = 0;				# quantile 0.5
			q5 = 1.644854;		# quantile 0.95
			q6 = 1.959964;		# quantile 0.975
			q7 = 2.575829;		# quantile 0.995
			quantiles = predict(rk.spline, c(q1, q2, q3, q4, q5, q6, q7))
			x_hat_low[i, "005"] = quantiles$y[1]
			x_hat_low[i, "025"] = quantiles$y[2]
			x_hat_low[i, "050"] = quantiles$y[3]
			x_hat_med[i] = quantiles$y[4]
			x_hat_upp[i, "950"] = quantiles$y[5]
			x_hat_upp[i, "975"] = quantiles$y[6]
			x_hat_upp[i, "995"] = quantiles$y[7]
		}
	} # for i
	
# Show the number of wrong cases (at most as many as time indices there are: N)
	cat("Number of cases with u*K'(u) - K(u) < 0:", wrong_count, ", smaller value:", wrong_min_value, "Summary:\n")
	print(summary(wrong_values))
	cat("Number of cases with K\"(u) < 0:", K_der2_neg_count, ", smaller value:", K_der2_min_value, "Summary: \n")
	print(summary(K_der2_neg_values))
	
#    list(x_hat_low=x_hat_low, x_hat_med=x_hat_med, x_hat_upp=x_hat_upp)
	list(x_hat_low=x_hat_low, x_hat_med=x_hat_med, x_hat_upp=x_hat_upp, u=u, r=r, Fu=Fu, K=K, K_der=K_der, K_der2=K_der2)
}
#################################### PROBABILITY BANDS ########################################


#****************************** FUNCTIONS FOR PROBABILITY BANDS *******************************
#**********************************************************************************************
