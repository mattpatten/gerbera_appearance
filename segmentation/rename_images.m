function rename_images(imgAbbr, sourceDir)

% Searches specified directory (full or relative) for any picture files and renames them as requested, followed by an index.
% For example, 03b8c701be95817e00ba0c9d5525dd15.jpg becomes pin01.png
%
% In addition, this script also saves a label file specifying the image names within the directory.
%
% If image directory isn't provided, assumes pictures are in current directory.
%
% Created by Matt Patten on 27 Nov 2018


%get image directory
[outputDir] = get_dir(imgAbbr,'img');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
%else
%    error(['Output folder: ' outputDir ' already exists.']);     
end

%Add file separator (/) if specified directory doesn't already include one
if ~strcmp(sourceDir(end),filesep)
    sourceDir = [sourceDir filesep];
end

%extract strings of the filenames of any picture type within this directory
image_filenames = get_image_filenames(sourceDir);

%load all images (done first to avoid over-writing any similarly named files later on)
%for indexed images (that rely on a separate colormap to save space), convert them to regular RGB images
for idx=1:length(image_filenames)
    
    %extract image file information
    warning off;
    full_filename = [sourceDir image_filenames{idx}];
    info = imfinfo(full_filename);
    warning on;
    
    if startsWith('tiff',info.Format)
        [flower(idx).img, map] = imread(full_filename,'Info',info);
    else
        try
            [flower(idx).img, map] = imread(full_filename);
        catch
            disp(['Error loading file: ' full_filename]);
        end
    end
    if ~isempty(map) 
        flower(idx).img = ind2rgb(flower(idx).img, map);
    end
end

%re-write image files with new names
for idx=1:length(image_filenames)
    labels{idx} = sprintf('%s%04i',imgAbbr,idx); %create new name
    imwrite(flower(idx).img, [outputDir labels{idx} '.png'], 'png'); %regardless of file type, saves as png file
end

save([outputDir 'labels_' imgAbbr '.mat'], 'labels');

end