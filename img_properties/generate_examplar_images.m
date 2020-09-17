function generate_examplar_images(imgAbbr)

% This loads the properties table and then ranks the scores on each property from lowest and highest
% and displays the gerberas that reflect these scores
%
% Input:
%   imgAbbr - The dataset abbreviation.
%
% Output:
%   A folder of images in output > imgAbbr > Exemplars with the top and bottom ranked flower images.
%
% Created by Matt Patten
% Created in Jan 2020.


%how many flowers do we want to display
numPlaces = 10;

dataDir = get_dir(imgAbbr,'data');

%load label and property details from labels file
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels');

for i=1:length(header)
    create_examplar(imgAbbr, labels, imgProperties(:,i), header{i}, numPlaces);
end

end