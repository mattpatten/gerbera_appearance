function [largestClass] = separateDistributions(vals, n)

% Takes values along a number line, and splits the data (via clustering) into n groups, returning the 
% indexes of the group that had the largest mean values.
%
% Inputs:
%    vals - the data to be separated into groups, a vector of real numbers.
%     n   - the number of groups to split the data into.
%
% Outputs:
%    largerClass - the indexes of values that belonged to the class with larger values
%
% Created by Matt Patten in Jan, 2019

% Run clustering to identify if flower is one colour or two
clusters = perform_kmeans(vals, n, 0); %perform k-means clustering

%get the mean value of each class
for i=1:n
    meanVal(i) = mean(vals(clusters==i));
end

%find the class with the highest mean
[~, maxClass] = max(meanVal);

%return the indexes of members of the class that contained larger values
largestClass = find(clusters==maxClass);

end