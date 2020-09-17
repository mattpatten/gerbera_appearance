function yhat = logistic2(params,x)

% Design matrix for logistic function.
% Returns vector of y-values corresponding to vector of
% x-values for the 2-element parameter vector "params". 
%
% CC 12.12.99

yhat = 1./(1+exp(-(x-params(1))./(params(2)/1.15)));

% init_params = [10 1]; % first guess to start iteration ...
% initial_fit = nlinfit(x,raw_probs,'logistic2',init_params)

%yhat is probably the estimated value of the slope given the threshold/slope parameters it is sent, at a particular value x
%we can then use yhat and compare that to the actual data point to see how close the estimate is

end
