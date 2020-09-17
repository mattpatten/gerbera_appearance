function [image_filenames] = get_image_filenames(sourceDir)

% For a specified directory, extracts the names of any image files located within.
%
% Created by Matt Patten
% Created in Feb 2020


image_types = {'tif','tiff','jpg','jpeg','png','bmp','gif'};

%if stimulus directory not specified, get current directory
if ~exist('sourceDir','var') 
    sourceDir = [pwd filesep];
end

%Add file separator (/) if specified directory doesn't already include one
if ~strcmp(sourceDir(end),filesep)
    sourceDir = [sourceDir filesep];
end

%get names of all images found within the directory
image_details = [];
for i=1:length(image_types)
    image_details = [image_details; dir(sprintf('%s*.%s',sourceDir,image_types{i}))];
end

% Extracts filename only
image_filenames = cell(1,length(image_details));
for idx=1:length(image_details)
    image_filenames{idx} = image_details(idx).name;
end

end