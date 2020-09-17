function [image] = applyMask(image,mask)

%Use the mask image to convert any regions outside the flower to NaN
%
% Inputs:
%    image - An image matrix, any number of dimensions.
%    mask  - A 2D logical matrix the same size, or the same size as the first two dimensions of the image.
%    
% Output:
%    image - The inputted image file where anywhere outside of the mask is changed from 0 to NaN 
%            across all provided dimensions.
%
% Created by Matt Patten, Oct 2018
% Modified, MP, Jan 2019: Removed need for conversion to double

if ~islogical(mask)
    mask = imbinarize(mask); 
end

if isequal(size(image),size(mask)) %if same dimensions of image and mask

    image(mask==0)=NaN; %replace anything outside the mask with NaNs

elseif size(image,1)==size(mask,1) && size(image,2)==size(mask,2) %xSize and ySize are consistent

    repMask = repmat(mask,1,1,size(image,3)); %replicate to the right number of dimensions
    image(repMask==0)=NaN; %replace anything outside the mask with NaNs

else
    error('Image and mask have different dimensions.');
end
    
end