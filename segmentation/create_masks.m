function create_masks(imgAbbr,pixSize)

% Uses a region growing algorithm to separate foreground/background between 
% the Gerbera and its surrounds.
%
% Inputs:
%    imgAbbr - Abbreviated description of the dataset of images to use.
%    pixSize - The length, in pixels, of each side of the image at the end of preprocessing
%
% Outputs:
%    A series of masks saved as image files: ones as foreground and zero as background.
%
% Created by Matt Patten in Dec, 2018

%parameters
tolerance = 20;

%file structure
[maskDir, imgDir] = get_dir(imgAbbr,'mask','img'); %directory for masks (located within image directory)
if ~exist(maskDir,'dir'), mkdir(maskDir); end %generate directory if not already existing

%load image details from labels file (or generate if not already there)
try load([imgDir 'labels_' imgAbbr '.mat'],'labels');
catch
    image_details = dir(sprintf('%s*.%s',imgDir,'png'));
    labels = cellfun(@(x) regexp(x,'.*(?=.png)','match','once'),{image_details.name},'UniformOutput',0);
    save([imgDir 'labels_' imgAbbr '.mat'], 'labels');
end

for flowerIdx = 1:length(labels)

    fprintf(['\nCreating mask for flower #' num2str(flowerIdx)]); %update user
    
    %load image and mask for this flower
    gerb.RGB  = imread([imgDir  labels{flowerIdx}  '.png']); %load image

    [xSize, ySize, ~] = size(gerb.RGB); %get image properties
    
    %Reduce image size to make processing time more manageable
    resized = false;
    while min(xSize,ySize) >= (2*pixSize)
        gerb.RGB = imresize(gerb.RGB,0.5);
        [xSize, ySize, ~] = size(gerb.RGB);
        resized = true;
    end
    if resized, imwrite(gerb.RGB,  [imgDir  labels{flowerIdx}  '.png'], 'png'); end %resave RGB image to new size
        
    %indexes of corners
    cnrs = [  1     1  ; ...
              1   ySize; ...
            xSize   1  ; ...
            xSize ySize];
    
    xCnrs = cnrs(:,1);
    yCnrs = cnrs(:,2);

    %{
    %% Old (significantly slower) algorithm for region growing to find flower edges.
    
    %convert to greyscale
    gerb.grey = rgb2gray(gerb.RGB);

    %define threshold (cutoff) for region growing algorithm (number of grey values) based on variance in pixels around the edge
    threshold_rg = decide_rg_threshold(gerb.grey);
    
    %% Perform region growing algorithm from each corner
    rg_all = zeros(xSize,ySize); %pre-allocate / reset
    for i=1:length(xCnrs)
        rg = regiongrowing(double(gerb.grey)/255, xCnrs(i), yCnrs(i), threshold_rg/255); %grow region from corner
        rg_all = or(rg_all,rg); %combine this with other corners
    end
    rg_all = ~rg_all; %anything that ~wasn't~ selected by the 4 region growing algorithm is part of the flower
    %}
    
    %% Magic wand method for region growing
    rg_all = zeros(xSize,ySize); %pre-allocate / reset
    for i=1:length(xCnrs)
        rg = magicwand2(gerb.RGB, tolerance, xCnrs(i), yCnrs(i)); %grow region from corner
        rg_all = or(rg_all,rg); %combine this with other corners
    end
    rg_all = ~rg_all; %anything that ~wasn't~ selected by the 4 region growing algorithm is part of the flower
    
    %select largest cluster and remove all else (borders, watermarks, etc)
    gerb.mask = get_biggest_cluster(rg_all);
    
    %alternative approach: binarize greyscale image, remove holes and then select largest cluster in image
    %naturalImg_mask = get_biggest_cluster(imfill(imbinarize(gerb.grey,'global'),'holes')); 

    %Save mask
    if ~isempty(labels{flowerIdx})
        imwrite(gerb.mask, [maskDir labels{flowerIdx} 'm.png'], 'png');
    end
end

end