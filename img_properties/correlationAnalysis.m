function [X, fitY, coeffs, r, p] = correlationAnalysis(Y)

% Examines 1D data and fits a line to it, returning the coeffecients (y = ax + b) and Pearson correlation values.
%
% Inputs:
%    Y      - A second column vector of the second variable to be correlated.
%
% Outputs:
%    X      - A vector of the x-values used in this analysis (after NaN trimming)
%    fitY   - The best fit of the linear function at each value of X.
%    coeffs - Two values indicating slope and intercept (y = ax + b), respectively.
%    r      - Pearson correlation of linear fit to data.
%    p      - Significance of Pearson correlation.
%
% Created by Matt Patten in Feb, 2019

%two input alternative (if given both x and y columns, rather than starting from 1 for x)
%corrData = [xvals yvals]; %put two columns - one as an incremental x value, the other as our variable of interest (e.g., hue)
%corrData = corrData(all(~isnan(corrData),2),:); %trim to relevant data (remove any row with a nan)

Y = Y(~isnan(Y)); %remove any NaN entries
X = (1:length(Y))'; %get x-values
corrData = [X Y]; %combine
coeffs = polyfit(corrData(:,1), corrData(:,2), 1); %fit linear function to data (y = ax + b)

fitY = polyval(coeffs, X); %get values at y for each x given the above coefficients
[tmp_r, p_r] = corr(corrData);
r = tmp_r(1,2); %get correlation, i.e., R^2 value
p = p_r(1,2);   %get p-value of correlation

end