function [cleanImg] = get_biggest_cluster(img)

% Goes through a boolean matrix and identifies the patch/cluster that is larger than all of the others,
% returning ones for that cluster and zeros for all the other smaller clusters.
%
% Inputs:
%    img      - A 2D binary image matrix.
%    
% Output:
%    cleanImg - The returned binary image where only the largest cluster is still present.
%
% Created by Matt Patten in Dec, 2018


%assigns labels to different clusters (if elements touch, it's considered part of the same cluster)
clusters = bwlabel(img); 

%replace non-cluster elements (the background / zeros) to nans so calculating the mode doesn't include this
clusters(clusters==0)=NaN; 

%find biggest cluster
cleanImg = clusters==mode(reshape(clusters,1,numel(clusters))); 

end