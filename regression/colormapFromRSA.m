function cols = colormapFromRSA

% this function provides a convenient colormap for visualizing
% dissimilarity matrices. it goes from blue to yellow and has grey for
% intermediate values.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

nCols = 256;
%% blue-cyan-gray-red-yellow with increasing V (BCGRYincV)
anchorCols=[0 0 1
            0 1 1
            .5 .5 .5 
            1 0 0
            1 1 0];

anchorCols_hsv=rgb2hsv(anchorCols);
incVweight=1;
anchorCols_hsv(:,3)=(1-incVweight)*anchorCols_hsv(:,3)+incVweight*linspace(0.5,1,size(anchorCols,1))';

anchorCols=hsv2rgb(anchorCols_hsv);

%% define color scale
nAnchors=size(anchorCols,1);
cols = interp1((1:nAnchors)',anchorCols,linspace(1,nAnchors,nCols));

end