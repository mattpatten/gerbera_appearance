function reorder_images(imgAbbr)

% Re-saves images from Holstein dataset in a more useful order (e.g., in order of participant preferences)
%
% Inputs:
%    imgAbbr - Abbreviated description of the dataset of images to use. 
%      Use: 'pref', 'pres' or 'orig' for order of preference or presentation, or original ordering.
%
% Output:
%    Image files and mat file for labels is saved in new image directory (under stimuli).
%
% Created by Matt Patten in 2018.


%get file structure
[rootDir, stimDir, imgDir, maskDir] = get_dir(imgAbbr,'root','stim','img','mask');
if ~exist(imgDir,'dir'), mkdir(imgDir); end
if ~exist(maskDir,'dir'), mkdir(maskDir); end

%load image names
if strcmpi(imgAbbr,'orig')
    %an overly complicated way to simply work out how many images there are to rename
    image_details = dir(sprintf('%s*.%s',[stimDir filesep 'used_images' filesep],'jpg')); 
    for i=1:length(image_details)
        labels_Holstein{i} = sprintf('H%02i',i-1); 
    end
    
else
    load(sprintf('%slabels_%s.mat',rootDir,imgAbbr)); %in order they were preferred/presented
    labels_Holstein = labels;
    clear labels;
end

%load images and masks
for flowerIdx = 1:length(labels_Holstein)
    % Imports image as such: gerbs (rows (x), cols (y), colour, image_idx)
    gerbs(:,:,:,flowerIdx)     = imread([stimDir 'used_images' filesep labels_Holstein{flowerIdx} '.jpg']);
    gerbsMask(:,:,:,flowerIdx) = imread([stimDir 'used_images' filesep 'masks' filesep labels_Holstein{flowerIdx} 'm.jpg']);
end

%save images and masks under new name/order
for flowerIdx = 1:length(labels_Holstein)
    imwrite(gerbs(:,:,:,flowerIdx)    ,sprintf('%s%s%02i.png',imgDir,imgAbbr,flowerIdx) ,'png');
    imwrite(gerbsMask(:,:,:,flowerIdx),sprintf('%s%s%02im.png',maskDir,imgAbbr,flowerIdx),'png');
    labels{flowerIdx} = sprintf('%s%02i',imgAbbr,flowerIdx); %create new name
end

%save names of new labels
save([imgDir 'labels_' imgAbbr '.mat'],'labels','labels_Holstein');

end