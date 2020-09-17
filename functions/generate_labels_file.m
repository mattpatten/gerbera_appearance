function generate_labels_file(imgAbbr, sourceDir)

% Creates labels file from images located within directory
%
% Created by Matt Patten
% Created in Feb 2020

%Add file separator (/) if specified directory doesn't already include one
if ~strcmp(sourceDir(end),filesep)
    sourceDir = [sourceDir filesep];
end

%extract strings of the filenames of any picture type within this directory
image_filenames = get_image_filenames(sourceDir);

%remove the extension (extract string elements before the full stop)
labels = cellfun(@(x) extractBefore(x,'.'), image_filenames, 'UniformOutput', 0);
%labels = cellfun(@(x) strcat(imgAbbr_new,extractBetween(x,imgAbbr_old,'.')), labels);

%save labels file
save([sourceDir 'labels_' imgAbbr '.mat'], 'labels');

end