function saveSegmentationFinal(imgAbbr)

% This function writes all images with the segmentation from the details stored in 
% the saved mat file. It's a good 'safety check' to make sure the segmentation is 
% still the way you think it is, and only takes a second to run.
%
% Inputs:
%   imgAbbr - The abbreviation of the dataset we're working with
%
% Output:
%   A directory filled with images with segmentation obtained from segmentation file
%
% Created by Matt Patten in January, 2019.


[imgDir, outputDir] = get_dir(imgAbbr,'img','output');
segDir = [outputDir 'segmentation_final' filesep];
if ~exist(segDir,'dir'), mkdir(segDir); end %if directory doesn't exist, create it

%load labels (filenames) and segmentation details
load([imgDir 'labels_' imgAbbr '.mat'],'labels');
try load([imgDir 'segmentation_manual.mat'],'segmentation'); %load manual segmentations (if file is there)
catch, load([imgDir 'segmentation.mat'],'segmentation');     %otherwise load automatically-generated version
end

for flowerIdx = 1:length(labels)

    %load image and mask for this flower
    gerb.RGB  = imread([imgDir labels{flowerIdx} '.png']); %load image
    filename = labels{flowerIdx};
    seg = segmentation(flowerIdx,:);
    
    %get properties
    [xSize, ySize, ~] = size(gerb.RGB);
    r = convert2polar([xSize ySize]);

    outputFile = [segDir filename '.png'];
    circles = or(r==seg(1), r==(seg(2)));
    imwrite(imoverlay(gerb.RGB,circles,'c'),outputFile,'png');

end
end