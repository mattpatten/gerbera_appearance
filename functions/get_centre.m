function [xCtr, yCtr] = get_centre(mask)

% Calculates centre of flower using median index in (logical) mask across both rows and columns
%
% Inputs:
%    mask - A binary image with ones inside object and zeros outside. 
%
% Outputs:
%    [xCtr, yCtr] - A variable pair containing the centre pixel for each image.
%
% Created by Matt Patten, 12/09/2018

%extract some image information
[xSize, ySize, ~] = size(mask); %get pixel dimensions of mask image
[imX, imY] = ndgrid(1:xSize,1:ySize); %generate index values for rows and columns

mask(mask==0)=NaN; %replace anything outside the mask with NaNs
xCtr = round(median(median((imX .* mask),2,'omitnan'),'omitnan')); %get median for each x-value (column) and again across all of these medians, rounding to nearest pixel value
yCtr = round(median(median((imY .* mask),1,'omitnan'),'omitnan')); %get median for each y-value (row) and again across all of these medians, rounding to nearest pixel value

end