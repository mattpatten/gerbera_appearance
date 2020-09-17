function remove_bad_images(imgAbbr)

% Remove any images where the algorithm has failed or were inappropriate (e.g., messy background)
%
% Inputs:
%    imgAbbr - Abbreviated description of the dataset of images to use.
%
% Outputs:
%    A series of flowers and masks moved to the subdirectory "rem".
%
% Created by Matt Patten in July, 2019


%parameters
threshold_mask = 2;        %percentage of pixels defined by the region growing algorithm before it should be deemed successful
colsToCheck = 5;           %number of columns (width, in pixels) to check around the border to identify if a natural image
similarity_threshold = 1;  %percentage of border pixels with the same grey value, below which assumes that masking procedure failed

%file structure
[maskDir, imgDir, remDir, remMaskDir] = get_dir(imgAbbr, 'mask', 'img', 'rem', 'remMask'); %directory for masks (located within image directory)

%load image details from labels file
load([imgDir 'labels_' imgAbbr '.mat'],'labels');
new_labels = labels; %updated label list as we're cutting images as we go (affecting the loop)

for flowerIdx = 1:length(labels)

    %load image and mask for this flower
    gerb.RGB  = imread([imgDir  labels{flowerIdx}  '.png']); %load image
    gerb.mask = imread([maskDir labels{flowerIdx} 'm.png']); %load mask
   
    [xSize, ySize, ~] = size(gerb.RGB); %get image properties
    
    %convert to greyscale
    gerb.grey = rgb2gray(gerb.RGB);
    
    % check if algorithm failed
    reggrow_fail = sum(sum(gerb.mask==0)) < (xSize*ySize)/100 * threshold_mask;
    if reggrow_fail, disp(['Removed flower #' labels{flowerIdx} ': Region growing failure.']); end
        
    %check background to see if natural image and get appropriate level depending on background consistency
    [natImg_flag, similarity] = check_for_natural_image(gerb.grey, colsToCheck, similarity_threshold); %if too many segments in first few columns, then likely natural image
    if natImg_flag, disp(['Removed flower ' labels{flowerIdx} ': Background similarity too low (' num2str(similarity) ').']); end
    
    %if natural background or algorithm failed, remove image
    if or(natImg_flag, reggrow_fail)
        if ~exist(remDir,'dir'),     mkdir(remDir);     end %generate directory for removed images if not already existing
        if ~exist(remMaskDir,'dir'), mkdir(remMaskDir); end %and also for their masks
        movefile([imgDir  labels{flowerIdx}  '.png'],[remDir     labels{flowerIdx}  '.png'],'f'); %move image
        movefile([maskDir labels{flowerIdx} 'm.png'],[remMaskDir labels{flowerIdx} 'm.png'],'f'); %move mask
        new_labels{flowerIdx} = ''; %remove from label
    end
end

%update label list (for images that have now been excluded)
labels = new_labels(~cellfun('isempty',new_labels));
save([imgDir 'labels_' imgAbbr '.mat'],'labels');

end