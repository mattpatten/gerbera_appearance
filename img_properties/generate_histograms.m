function generate_histograms(imgAbbr)

% This loads the properties table and then summarises each proponent of the table the data through histograms
% 
%
% Input:
%   imgAbbr - The dataset abbreviation.
%
% Output:
%   A folder of images in output > imgAbbr > Histograms with a histogram of each property.
%
% Created by Matt Patten
% Created in Feb 2020.


%how many flowers do we want to display
[dataDir, outputDir] = get_dir(imgAbbr,'data','output');
saveDir = [outputDir 'Histograms' filesep]; %define directory
if ~exist(saveDir,'dir'), mkdir(saveDir); end %if directory doesn't exist, create it

%load label and property details from labels file
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels');

for i=1:length(header)
    f = figure;
    histogram(imgProperties(:,i));
    xlabel('Property Value');
    ylabel('Freq');
    title(header{i},'Interpreter','None');
    box off;
    saveas(f,[saveDir header{i} '.png']);
    saveas(f,[saveDir header{i} '.fig']);
    close;
end

end