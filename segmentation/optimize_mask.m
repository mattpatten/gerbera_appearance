function [xCtr, yCtr, rad] = optimize_mask(mask, saveDir, label)

% This function uses fminsearch to identify the flower face, that is, the centre 
% and radius of a circle which contains the most amount of pixels in the mask.
%
% Inputs:
%   mask    - A logical matrix describing parts of the flower (ones) and parts outside the flower (zeros)
%   saveDir - Directory for saving output. 
%   label   - Flower name/index for save file.
%
% Outputs:
%   xCtr    - The optimized x-position of flower centre.
%   yCtr    - The optimized y-position of flower centre.
%   rad     - The optimal radius of the flower face.
%
% Created by Matt Patten in Dec 2018


%get properties
[xSize, ySize] = size(mask);

%get initial parameters (centre)
[init_xCtr, init_yCtr] = get_centre(double(mask));

%get initial parameters (radius)
for i = 1:100 %nsteps
    radsteps = linspace(1,min([xSize ySize])/2);
    radii(i) = calcMaskPixels(log([init_xCtr, init_yCtr radsteps(i)]), mask); %go through all radii and select one that is best for this (inital) centre
end
[~, init_rad_idx] = min(radii); %get radius size corresponding to the local minimum
init_rad = radsteps(init_rad_idx);
disp('Radius initialized');

%save co-ordinates to later display the initial parameters (and see how much it has changed from the algorithm)
r = convert2polar([xSize ySize], [init_xCtr init_yCtr]);
init_ctr = or(r==0,r==1);
init_circ = r==round(init_rad);

%find centre and radius of circle that contains maximum number of pixels in the mask
init_params = log([init_xCtr, init_yCtr, init_rad]); %log as algorithm works better with values at a smaller scale
param = fminsearch(@(param) calcMaskPixels(param, mask), init_params);

%save final variables
xCtr = round(exp(param(1))); %undo log, round to nearest pixel/radius
yCtr = round(exp(param(2)));
rad  = round(exp(param(3)));

%generate logical matrices showing centre and radius locations
r = convert2polar([xSize ySize], [xCtr yCtr]);
ctr = or(r==0,r==1);
circ = r==rad;

%display results
imwrite(imoverlay(imoverlay(imoverlay(imoverlay(mask,init_ctr,[0.8 0.8 0.8]),init_circ,[0.8 0.8 0.8]),ctr,'r'),circ,'r'),[saveDir 'ctr_' label '.png'],'png');

end


function [invNumMaskPix] = calcMaskPixels(p, mask)

    % parameters: xCtr, yCtr, radius

    %get matrix of the distance of all pixels from these centre co-ordinates
    r = convert2polar(size(mask), round(exp([p(1) p(2)]))); 
    
    %get number of pixels within and including this radius
    circle = r<=round(exp(p(3))); 
    
    %find number of values in this radius that are also in the mask and subtract pixels that aren't. Take recipricol (as we're finding a minimum).
    invNumMaskPix = 1/(nnz(and(circle,mask)) - nnz(and(circle,~mask))); 
    
    %a negative number indicates more pixels outside of the mask than inside the mask and should be excluded (esp. since we're inverting and finding minimum)
    if invNumMaskPix<=0, invNumMaskPix=NaN; end 
end
