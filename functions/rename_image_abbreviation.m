function rename_image_abbreviation(oldImgAbbr, newImgAbbr, sourceDir)

% Takes a stimulus set and copies all files to a new folder with the new imgAbbr, updates the labels and filenames
%
% Created by Matt Patten
% Created in Feb 2020

outputDir = get_dir(newImgAbbr,'img');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

%Add file separator (/) if specified directory doesn't already include one
if ~strcmp(sourceDir(end),filesep)
    sourceDir = [sourceDir filesep];
end

%extract strings of the filenames of any picture type within this directory
%image_filenames = get_image_filenames(sourceDir);

%load and rename labels
load([sourceDir 'labels_' oldImgAbbr '.mat'],'labels');
old_labels = labels;
labels = cellfun(@(x) strcat(newImgAbbr,extractAfter(x,oldImgAbbr)), labels, 'UniformOutput', 0);
save([outputDir 'labels_' newImgAbbr '.mat'], 'labels');

%copy and rename files
for i=1:length(old_labels)
    copyfile([sourceDir old_labels{i} '.png'],[outputDir labels{i} '.png']);
end

%copy segmentation details to new directory
try copyfile([sourceDir 'segmentation_manual.mat'],[outputDir 'segmentation_manual.mat']); 
catch, copyfile([sourceDir 'segmentation.mat'],[outputDir 'segmentation.mat']); 
end

end