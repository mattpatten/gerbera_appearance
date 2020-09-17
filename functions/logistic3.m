function yhat = logistic3(params,x)

% Design matrix for logistic function.
% Returns vector of y-values corresponding to vector of
% x-values for the 3-element parameter vector "params"
% where the 3rd parameter is the miss rate.
%
% CC 28.09.01

yhat = params(3)/2 + (1-params(3))./(1+exp(-(x-params(1))./(params(2)/1.15)));

end
