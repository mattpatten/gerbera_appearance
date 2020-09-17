function generate_segmentation_subset(imgAbbr, sourceDir)

% Opens the original segmentation file and its corresponding labels file, and compares with a reduced labels file for 
% a subset of the stimuli, to similarly reduce the segmentation file to just this new subset of stimuli.
%
% Created by Matt Patten
% Created in Feb 2020

%Add file separator (/) if specified directory doesn't already include one
if ~strcmp(sourceDir(end),filesep)
    sourceDir = [sourceDir filesep];
end

%load label details from full dataset
[imgDir] = get_dir(imgAbbr,'img');
load([imgDir 'labels_' imgAbbr '.mat'],'labels');
labels_full = labels;
clear labels;

%load label details from data subset
load([sourceDir 'labels_' imgAbbr '.mat'],'labels');
labels_subset = labels;
clear labels;

%load segmentation details
try 
    load([imgDir 'segmentation_manual.mat'],'segmentation'); %load manual segmentations (if file is there)
    disp('Loading most recent segmentation file with some manual corrections.');
    manual_text = '_manual';
catch
    load([imgDir 'segmentation.mat'],'segmentation');     %otherwise load automatically-generated version
    disp('Loading original (auto) segmentation file.');
    manual_text = '';
end

%find subset of label images and update segmentation accordingly
[~, idx] = ismember(labels_subset,labels_full);
segmentation = segmentation(idx,:);

%save in image subset directory
save([sourceDir 'segmentation' manual_text '.mat'],'segmentation'); %Careful - this will overwrite existing files

end