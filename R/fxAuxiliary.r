#**********************************************************************************************
#*********************************** AUXILIARY FUNCTIONS **************************************

# INDEX
# =====
# bound


# Function that bounds an array to a given bound.
bound = function(x, bound, type)
## Input parameters:
## x:       Value to be bounded
## bound:   Value of the bound
## type:    Type of bound: MIN (upper bound), MAX (lower bound), BOTH (both lower and upper bound.
##          In this case, the bound is assumed to be positive)
{
	if (toupper(type) == "MIN")
	{
		ind = which(x > bound);
		if (length(ind) > 0)
		{
			cat(paste("---OVERFLOW (1-POS, bound=", bound, ")---\n", sep=""))
			cat("length(ind)= ", length(ind), "\n")
			cat("Summary of values before bounding\n")
			print(summary(x[ind]))
		}
		x[ind] = bound;
	}
	else if (toupper(type) == "MAX")
	{
		ind = which(x < bound);
		if (length(ind) > 0)
		{
			cat(paste("---OVERFLOW (2-NEG, bound=", bound, ")---\n", sep=""))
			cat("length(ind)= ", length(ind), "\n")
			cat("Summary of values before bounding\n")
			print(summary(x[ind]))
		}
		x[ind] = bound;
	}
	else
	{
		# bound is assumed to be positive.
		indlow = which(x < -bound);
		indupp = which(x > bound);
		if (length(indlow) > 0)
		{
			cat(paste("---OVERFLOW (3-NEG, bound=", bound, ")---\n", sep=""))
			cat("Summary of values before bounding\n")
			print(summary(x[indlow]))
		}
		if (length(indupp) > 0)
		{
			cat("---OVERFLOW (3-POS)---\n")
			cat("Summary of values before bounding\n")
			print(summary(x[indupp]))
		}
		x[indlow] = -bound;
		x[indupp] = bound;
	}
	return(x)
}


#*********************************** AUXILIARY FUNCTIONS **************************************
#**********************************************************************************************
