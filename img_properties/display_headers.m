function display_headers(imgAbbr)

% This displays on screen each header from the properties table with numbering.
%
% Created by Matt Patten
% Created in Feb 2020


%load label and property details from labels file
dataDir = get_dir(imgAbbr,'data');

%load image properties table
load([dataDir 'properties_table_' imgAbbr '.mat'],'header');

%display image properties with indexes
disp(' ');
disp([num2cell((1:length(header))') header']);

end