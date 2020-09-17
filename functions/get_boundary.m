function [border] = get_boundary(img,d)

% This identifies and returns the regions within a set distance from the outside edge of the flower.
% It does this by grabbing a box around each pixel and seeing if there's a NaN in that box 
% (i.e., a region already masked and identified as a non-flower pixel). If there is, it's close
% to the edge of the flower.
% 
% Inputs:
%    img    - A 2D image matrix, already masked with NaN's for pixels outside of the mask.
%    d      - The distance, in pixels, from the border of the flower to examine.
%
% Outputs:
%    border - A 2-D logical matrix of size 'img' that identify pixels that are within 
%             a set distance to the border.
%
% Created by Matt Patten in Nov 2018

%Initialize / get properties
[xSize, ySize] = size(img);
border = zeros(xSize,ySize);
boxWidth = 2*d+1; %how many pixels wide we will be checking for the neighbourhood

%get neighbourhood values of each pixel
boxCols = im2col(img,[boxWidth boxWidth]); 

%find any elements that have a NaN in the neighbourhood
nanEntries = any(isnan(boxCols),1); 

%put back into original shape
border((1+d):(end-d),(1+d):(end-d)) = reshape(nanEntries, xSize-2*d,ySize-2*d); 

%insert mask, convert to logical
border(isnan(img))=0; 
border = logical(border);

end