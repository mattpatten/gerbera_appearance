function [gerbs, labels] = load_images(imgAbbr, isStruct)

% Loads gerbera images and their corresponding labels.
%
% Inputs:
%    imgAbbr  - Abbreviated description of the dataset of images to use.
%    isStruct - Whether the gerberas should be returned as a struct (one per image) or a 4-D array.
%
% Outputs:
%    gerbs    - The image matrices.
%    labels   - A cell array with the labels for each of the images.
%
% Created by Matt Patten in Dec 2018

[imgDir] = get_dir(imgAbbr,'img');

%load image labels
load([imgDir 'labels_' imgAbbr '.mat']);

for flowerIdx = 1:length(labels)
    
    if isStruct
        % Imports image as a struct, in the format:
        % gerbs(image index).RGB(x pos (col), y pos (row), colour). Note oddity that x=col and y=row.
        gerbs(flowerIdx).RGB = imread([imgDir labels{flowerIdx} '.png']);
        
    else
        % Imports image as a 4-D array, in the format:
        % gerbs(x pos (col), y pos (row), colour, image index).
        gerbs(:,:,:,flowerIdx) = imread([imgDir labels{flowerIdx} '.png']);
        
    end
end

end