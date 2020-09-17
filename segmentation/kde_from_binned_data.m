function [kde, p] = kde_from_binned_data(bins, binned_data, mult, alpha)

% Computes the kernel density estimate for already-binned data. 
% 'Binned data' refers to lists of values as follows:
%   1     7
%   2   210
%   3   311
%   .    .
%   .    .
%   .    .
%  150   42
%  whereas the kde function requires input such as [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, ....... , 150, 150, 150, 150, 150]
% 
% Inputs: 
%    bins        - The range of values that the binned data represent. e.g., 1:150 for 150 bins from 1 to 150.
%    binned_data - A vector containing the frequencies or proportions of each value (e.g., 1: 0.47, 2: 0.21, 3: 0.91 or alternatively 1: 23, 2: 91, 3: 44)
%    mult        - If the frequencies are proportions, a multiplier value can be used to convert each entry to a positive integer [default=1]
%    alpha       - The value of statistical significance for the area underneath the kernel smoothed function [default=0.01]
%
% Output:
%    kde         - A kernel density estimate for each bin value, scaled to the original data 
%    p           - Head and tail percentiles from the data, once unbinned (as specified by alpha). 
%
% Created by Matt Patten in Nov, 2018


%define default values for those that can be assumed
if ~exist('mult'),  mult = 1;     end
if ~exist('alpha'), alpha = 0.01; end

%remove data from bins data
unbinned_data = [];
for i=bins
    unbinned_data = [unbinned_data, repmat(i,1,round(binned_data(i) * mult))];
end

%compute kernel density estimate
try
    bw = 1.06 * nanstd(unbinned_data) * length(unbinned_data) ^ (-1/5); %Silverman's rule of thumb bandwidth estimator
    kde = ksdensity(unbinned_data, bins, 'Bandwidth', bw);
    kde = kde / max(kde) * max(binned_data); %scale kde pdf to original data
catch %if data is empty and kernel density estimation fails
    kde = zeros(1,length(bins));
end

%calculate percentiles
p = round(prctile(unbinned_data,[alpha,1-alpha]*100));

end