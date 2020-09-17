function sumsqerr = compute_fit_residuals(params,x,prop)

%for each of the x-values, get the value from the fit so that we can subsequently judge how close it is to our data
fit_prop = logistic4(params,x);

%the sum of squares of the statistical errors - specifically, the fit of the logistic function to our data, squared, and then summed for all data points
%this is accumulated for each data point across all runs, and is the output function of the fminsearch - i.e., this is what we are trying to find the minimum value of
sumsqerr = sum((prop-fit_prop).^2); 

end