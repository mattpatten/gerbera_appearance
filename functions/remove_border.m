function [newImg] = remove_border(img, border, mask)

% Sets both the pixels outside the flower (masked) and those around the edge (border) to zero from the 
% original image, returning the image borderless and without surround.
%
% Inputs:
%    img    - An image matrix.
%    border - A 2D logical matrix of size 'img', where ones are pixels close to the edge of the image object.
%    mask   - A 2D logical matrix of size 'img', where ones are inside the image object and NaN otherwise.
%    
% Output:
%    newImg - The original image, with pixels close to the border and outside the mask set to zero.
%
% Created by Matt Patten in Dec, 2018


img(border)=0; %remove pixels from the peripheral border region
img(isnan(mask))=0; %remove pixels from outside of the mask
newImg = img;

end