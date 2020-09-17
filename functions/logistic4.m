function yhat = logistic4(params,x)

% Design matrix for logistic function.
% Returns vector of y-values corresponding to vector of
% x-values for the 4-element parameter vector "params"
% where the 3rd and 4th parameters are the min/max values.
%
% CC  28.09.01
% MLP 05.02.19 - changed from 3-param (miss rate) to 4-param (min/max)

yhat = params(4) + (params(3)-params(4))./(1+exp(-(x-params(1))./(params(2)/1.15)));

end
