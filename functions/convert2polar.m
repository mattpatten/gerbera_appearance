function [r, theta] = convert2polar(imgSize, ctr)

% This returns a grid with the distance and angle of each pixel from the centre.
% If no centre location is given, it is assumed to be at the image centre (rather than defined centre)
% 
% Inputs:
%    imgSize - Two positive integers, in pixels, representing the x and y size co-ordinates
%    ctr     - Two positive integers describing the centre pixel of the image (in order x,y --> col,row) 
%
% Outputs:
%    r     - in pixels, the distance of each pixel to the centre
%    theta - the angle of the pixel from the centre
%
% Created by Matt Patten in Nov, 2018

%set default centre if not provided
if ~exist('ctr','var')
    ctr = [imgSize(1)/2 imgSize(2)/2]; 
end

%create grid
[X,Y] = ndgrid(linspace(-imgSize(1)/2,imgSize(1)/2,imgSize(1)),...
               linspace(-imgSize(2)/2,imgSize(2)/2,imgSize(2))); 

%shift centre if necessary
X = X-(ctr(1)-imgSize(1)/2); 
Y = Y-(ctr(2)-imgSize(2)/2);

%converts from cartesian to polar coordinates for all points within the specified grid
[theta,r] = cart2pol(X,Y); 

r = round(r); %give whole pixels
%theta = flipud(theta); %flip bottom half so orientations move anti-clockwise

end