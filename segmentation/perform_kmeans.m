function [cluster_labels, p, centroids, cent_dist] = perform_kmeans(data, num_centroids, noisy_pixels)

% Performs k-means clustering on the dataset.
% 
% Inputs:
%    data           - Data that is to be split into clusters, with each row being a single entry and each column 
%                     being a dimension for the purposes of clustering. That is, 3 columns is a point in 3D space.
%    num_centroids  - The number of clusters you want the data separated into.
%    noisy_pixels   - The number of initial pixels to exclude due to unreliable signal (e.g., too close to centre)
%
% Output:
%    cluster_labels - returns the cluster associated with each value in data, either as a class label (better 
%                     for logic/understanding) or specified as the centroid value (better for plots).
%    p              - The boundary values of the fit. Specifically, the first value of the second cluster to appear
%                     and the last value of the second-last cluster in the dataset.
%
% Created by Matt Patten in Nov, 2018


%perform k-means clustering
[cluster_labels, centroids, ~, cent_dist] = kmeans(data((noisy_pixels+1):end,:),num_centroids);
cluster_labels = padarray(cluster_labels, noisy_pixels, cluster_labels(1), 'pre'); %replace innermost (noisy) eccentricities with whatever cluster was identified first

%sort cluster classes (so 1 is closest to centre and n is farthest)
[~, order] = sort(sqrt(sum(centroids(:,1).^2,2))); %Pythagoras, i.e., Euclidian distance from zero.
tmp = NaN(size(data,1),1);
for km = 1:num_centroids
    tmp(cluster_labels==order(km))=km; %to keep as class labels
    %tmp(cluster_labels==order(km))=centroids(order(km)); %to keep as centroid value
end
cluster_labels = tmp;

%Get boundary values
p = [find(diff(cluster_labels)~=0, 1, 'first')+1 ...
     find(diff(cluster_labels)~=0, 1, 'last')]; %no +1 because I want the last value before it changes (+1 as diff starts on 2nd index, -1 to get value before it changes)

end