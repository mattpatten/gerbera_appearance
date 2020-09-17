function resize_images(imgAbbr, pixSize)

% This function uses fminsearch to identify the flower face, that is, the centre 
% and radius of a circle which contains the most amount of pixels in the mask.
%
% Inputs:
%    imgAbbr - Prefix to name each image file, and partial name of stored directory (images_imgAbbr)
%    pixSize - The length, in pixels, of each side of the image at the end of preprocessing
%
% Outputs:
%   The stimulus and mask images, all resized to a 1:1 aspect ratio (i.e., a square) with white
%   background and to a specific pixel size.
%
% Created by Matt Patten in Dec 2018


%file structure
[imgDir, maskDir, procDir, remDir, remMaskDir] = get_dir(imgAbbr,'img','mask','proc','rem','remMask');
if ~exist(procDir,'dir'), mkdir(procDir); end %generate directory if not already existing

%load image details from labels file
load([imgDir 'labels_' imgAbbr '.mat'],'labels');

for flowerIdx = 1:length(labels)

    %load image and mask for this flower
    gerb.RGB  = imread([imgDir  labels{flowerIdx}  '.png']); %load image
    gerb.mask = imread([maskDir labels{flowerIdx} 'm.png']); %load mask

    %update user
    disp(['Resizing flower #' num2str(flowerIdx)]);
    
    %get ideal properties of flower face
    [xCtr, yCtr, rad] = optimize_mask(gerb.mask, procDir, labels{flowerIdx});

    %if flower head is too small in the image 
    if (rad*2)<pixSize/1.65 %flower head takes up [1=100%, 1.5=75%, 1.75=62.5%, 2=50%] of final image size

        if ~exist(remDir,'dir'),      mkdir(remDir);     end %generate directory if not already existing
        if ~exist(remMaskDir,'dir'),  mkdir(remMaskDir); end %and for their mask
        movefile([imgDir  labels{flowerIdx}  '.png'],[remDir     labels{flowerIdx}  '.png'],'f'); %move image
        movefile([maskDir labels{flowerIdx} 'm.png'],[remMaskDir labels{flowerIdx} 'm.png'],'f'); %move mask
        labels{flowerIdx} = '';
        disp([' REMOVED! Flower too small in image. (Diam: ' num2str(round(rad*2)) ', Image size: ' num2str(pixSize) ')']);

    else
        %fix aspect ratio - make image square (pad and/or crop)
        gerb.mask = make_image_square(gerb.mask, xCtr, yCtr, rad);
        gerb.RGB  = make_image_square(gerb.RGB,  xCtr, yCtr, rad);
        
        %get new stimulus size
        [xSize, ySize, nColours] = size(gerb.RGB);
        
        if xSize~=ySize
            disp(['WARNING: Image dimensions are not consistent: ' num2str(xSize) 'x' num2str(ySize)]);
        end
        
        %remove anything outside the mask from the image (standardizes background)
        gerb.RGB(~gerb.mask(:,:,ones(1,nColours))) = 255;
        
        %rescale image to specified size
        gerb.mask = imresize(gerb.mask,pixSize/xSize,'OutputSize',[pixSize pixSize]);
        gerb.RGB  = imresize(gerb.RGB, pixSize/xSize,'OutputSize',[pixSize pixSize]);
        
        %Save re-sized images
        imwrite(gerb.RGB,  [imgDir  labels{flowerIdx}  '.png'], 'png'); %save RGB image
        imwrite(gerb.mask, [maskDir labels{flowerIdx} 'm.png'], 'png'); %save mask
    end
end

%update label list (for images that have now been excluded)
labels = labels(~cellfun('isempty',labels));
save([imgDir 'labels_' imgAbbr '.mat'],'labels');

end
