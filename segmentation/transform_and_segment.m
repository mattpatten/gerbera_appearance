function transform_and_segment(imgAbbr)

% This function generates the necessary transformations to the Gerbera images and then uses these to find 
% the boundary between (1) the disk and trans florets, and (2) the trans and ray florets.
% 
% Inputs:
%    imgAbbr - Abbreviated description of the dataset of images to use.
%    
% Output:
%    The data is saved to a segmentation file within the directory where the images are stored. Data is a 
%    2-column vector storing the radial distance from the centre to these boundaries.
%
% Created by Matt Patten, Dec 2018

warning off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          File I/O           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up file I/O, load images & masks
[procStimDir, imgDir, maskDir, outputDir] = get_dir(imgAbbr,'procStim','img','mask','output');

%load image details from labels file
load([imgDir 'labels_' imgAbbr '.mat'],'labels');

for flowerIdx = 1:length(labels)

    disp(['Performing transformations and segmentation on flower #' num2str(flowerIdx)]); %update user
    
    %load image and mask for this flower
    gerb.RGB  = imread([imgDir  labels{flowerIdx}  '.png']); %load image
    gerb.mask = double(imread([maskDir labels{flowerIdx} 'm.png'])); %load mask

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Perform transformations   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Parameters
    boundary_width = 5; %in pixels, number of peripheral pixels of flower edge to consider as border

    % Convert from RGB to CIE L*a*b* and greyscale colourspaces
    disp('Converting to different colourspaces...'); %update status

    gerb.LAB  = rgb2lab(gerb.RGB);
    gerb.grey = rgb2gray(gerb.RGB); %conversion removes hue and saturation, keeps luminance - applies weighted average: 0.29R + 0.59G + 0.11B.

    % Apply mask to Gerberas
    disp('Applying mask....');

    %Use the mask image to convert any regions outside the flower to NaN
    gerb.greyMasked = applyMask(double(gerb.grey), gerb.mask); %greyscale
    for cc=1:size(gerb.RGB,3)
        gerb.LABMasked(:,:,cc) = applyMask(gerb.LAB(:,:,cc), gerb.mask); %coloured
    end

    %% Get flower boundary
    disp('Finding boundary pixels...');

    gerb.border = get_boundary(gerb.greyMasked, boundary_width);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Phase congruency     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Computing phase congruency...'); %update user

    %define parameters
    p.pc.k = 15; %checked between 5-20, 10-20 is the sweet spot
    p.pc.nscale = 2; %checked between 2-6: 2 is the best, 3 isn't so bad either (int only)
    p.pc.minWavelength = 4; %best to leave centre but get rid of rest

    %perform phase congruency analysis
    img_edge = phasecong3(gerb.grey, 'k', p.pc.k, 'nscale',p.pc.nscale, 'minWaveLength', p.pc.minWavelength);
    gerb.phaseCong = remove_border(img_edge, gerb.border, gerb.greyMasked);  %remove outer edge of flower (edges)

    clear img_edge;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Radial Segmentation    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %define the set of images we are going to use for this analysis
    img.rs.analysis = gerb.LABMasked;
    img.rs.display  = gerb.RGB;
    img.rs.phase    = gerb.phaseCong; %phase congruency from previous analysis

    %median
    p.rs.num_pixels = 1; %the width of the annulus of each bin
    p.rs.remove_pixels = 3; %remove the innermost number of circles due to lack of data
    p.rs.multiplier = 1000; %to turn proportions into frequencies for kde
    p.rs.sigValue = 0.05; %endpoint values to consider, each tail

    %k-means
    p.rs.num_centroids = 4; %number of centroids in k-means clustering algorithm

    %output
    p.rs.dispFig = 1; %to show output figures or not

    %generate file structure
    names.rs.image = labels{flowerIdx};
    names.rs.analysis = 'segmentation_auto';
    names.rs.imgDir = [procStimDir 'images_' imgAbbr filesep];
    names.rs.dir = [outputDir names.rs.analysis filesep];
    if ~exist(names.rs.dir,'dir'), mkdir(names.rs.dir); end %if directory doesn't exist, create it

    %perform analysis
    segmentation(flowerIdx,:) = radial_segmentation(img.rs, p.rs, names.rs);

end

%Save segmentation in image directory
save([imgDir 'segmentation.mat'],'segmentation');

end