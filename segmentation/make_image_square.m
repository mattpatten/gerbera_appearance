function [img] = make_image_square(img, xCtr, yCtr, rad)

% This function uses the centre coordinates and the radius of the flower face, then generates a 
% buffer/border around this and then crops/pads the image to convert the image into a 1:1 aspect 
% ratio (a square).
% 
% Inputs:
%   img     - The image to be modified into a 1:1 aspect ratio (i.e., a square)
%   xCtr    - The optimized x-position of flower centre.
%   yCtr    - The optimized y-position of flower centre.
%   rad     - The optimal radius of the flower face.
%
% Outputs:
%   img     - The input image, now modified into a 1:1 aspect ratio (i.e., a square)
%
% Created by Matt Patten in Dec 2018


rad_padding = 0.3; % the proportion of padding we want around optimal radius to ensure edge of petals are included and for minor border

padValue = 'replicate'; %the value (numerical value representing colour) or method (circular/replicate/symmetric) to pad array with

rad = rad + round(rad * rad_padding); %increase radius to ensure it includes edges of petals

[xSize, ySize, ~] = size(img); %get image properties

%compute extremities of new circle
xMin = xCtr - rad;
xMax = xCtr + rad;
yMin = yCtr - rad;
yMax = yCtr + rad;

if xMax < xSize %do max first so cropping doesn't affect indexes
    img = img(1:xMax,:,:); %crop any excess values (beyond radius)
else
    img = padarray(img,[(xMax-xSize) 0 0],padValue,'post'); %extend image with zeros
end

if yMax < ySize
    img = img(:,1:yMax,:); %crop any excess values (beyond radius)
else
    img = padarray(img, [0 (yMax-ySize) 0],padValue,'post'); %extend image with zeros
end

if xMin > 0
    img = img(xMin:end,:,:); %crop any excess values (beyond radius)
else    
    img = padarray(img, [abs(xMin)+1 0 0],padValue,'pre'); %extend image with zeros (+1 to include pixel for zero)
end

if yMin > 0
    img = img(:,yMin:end,:); %crop any excess values (beyond radius)
else
    img = padarray(img, [0 abs(yMin)+1 0],padValue,'pre'); %extend image with zeros (+1 to include pixel for zero)
end

end