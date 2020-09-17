function [maskSeg1, maskSeg2, kmeans_idxs] = segmentFloretsRadially(radialMean, mask)

% This uses L*a*b* colours to segment a floret region into two parts, based on a radial analysis.
%
% Inputs:
%    radialMean - A cell array containing 3xN matrices (LAB/RGB/HSV by number of eccentricities), indicating the mean colour for each channel at each eccentricity.
%    mask       - A logical matrix separating pixels inside and outside of the flower.
%
% Outputs:
%    maskSeg1   - Splits a region of the flower into two components based on clustering of radial (LAB) values.
%
% Created by Matt Patten in Jan, 2019


warning off;

%initial properties
numCentroids = 2;
edge_pixels_to_remove = 0; %not fully functioning within segments
[xSize, ySize] = size(mask);
r = convert2polar([xSize ySize]); 

%cluster analysis
kmeans_classes = perform_kmeans(radialMean, numCentroids, edge_pixels_to_remove); %perform k-means clustering
clusterChange = [NaN; (diff(kmeans_classes))]; %when there's a change in cluster label
kmeans_idxs = find(clusterChange==1 | clusterChange==-1)'; %get location of those changes

%find all eccentricities for one class/segment, get its circle and combine with mask to ensure it doesn't go outside flower
maskSeg1 = and(arrayfun(@(x) any(x==find(kmeans_classes==1)),r),mask);
maskSeg2 = and(arrayfun(@(x) any(x==find(kmeans_classes==2)),r),mask);
        
warning on;
end