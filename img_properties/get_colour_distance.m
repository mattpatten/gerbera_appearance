function [colourDist] = get_colour_distance(img, mask_seg1, mask_seg2)

% Computes the mean colour in each of two segments of the same image, and calculates and returns the Euclidean
% distance between the two colours.
%
% Inputs:
%    img       - An image to be processed and have its colour compared across two different segments (of size N x M x 3)
%    mask_seg1 - A logical matrix of size N x M x 1 operating as a first mask for its average colour to be computed.
%    mask_seg2 - Same, for a second mask to be compared against the mean colour of the first mask.
%
% Outputs:
%    colourDist - The Euclidean distance between the colours of the two segments.
%
% Created by Matt Patten in Feb, 2019


%get mean L*/a*/b* values for each of the segments
imgMean_seg1 = [nanmean(applyMask(img(:,:,1), mask_seg1), 'all') ...
                nanmean(applyMask(img(:,:,2), mask_seg1), 'all') ...
                nanmean(applyMask(img(:,:,3), mask_seg1), 'all')];

imgMean_seg2 = [nanmean(applyMask(img(:,:,1), mask_seg2), 'all') ...
                nanmean(applyMask(img(:,:,2), mask_seg2), 'all') ...
                nanmean(applyMask(img(:,:,3), mask_seg2), 'all')];

%get Euclidean distance between the mean values of the two segments
colourDist = pdist([imgMean_seg1; imgMean_seg2]);

end
