function create_examplar(imgAbbr,labels,score,statLabel,numPlaces)

% This takes a statistic, sorts it and displays a figure showing the flowers that have the highest
% and lowest values on this statistic.
%
% Inputs: (NB: score and labels must be vectors of the same lengths)
%   imgAbbr      - Dataset abbreviation label.
%   score        - The statistic to be ranked and displayed.
%   labels       - The labels for the images to be displayed.
%   statLabel    - The name of the statistic associated with this score.
%   numPlaces    - The number of images we want to display on each end (e.g., 10 highest and 10 lowest ranked flowers)
%
% Created by Matt Patten
% Created on 28/9/2018


% File I/O
[outputDir, imgDir] = get_dir(imgAbbr,'output','img');

saveDir = [outputDir 'Exemplars' filesep]; %define directory
if ~exist(saveDir,'dir'), mkdir(saveDir); end %if directory doesn't exist, create it

% Initialization
placeLabels = {'','2nd','3rd','4th','5th','6th','7th','8th','9th','10th'}; %exclude 1st as "1st Lowest" is nonsensical

%sort score
[val, idx] = sort(score); 

%{
%optional - remove zeros from display (i.e., if looking at gradient change, remove values that don't have a gradient at all)
if min(val)==0
    idxs_to_keep = val~=0;
    val = val(idxs_to_keep);
    idx = idx(idxs_to_keep);
end
%}
numFlowers = length(val);

%get values and indexes
lowestValues  = val(1:numPlaces);
lowestIdxs    = idx(1:numPlaces);
highestValues = flip(val((numFlowers-numPlaces):numFlowers));
highestIdxs   = flip(idx((numFlowers-numPlaces):numFlowers));

middle    = round(numFlowers/2-numPlaces/2);
midValues = val(middle:(middle+numPlaces-1));
midIdxs   = idx(middle:(middle+numPlaces-1));

f = figure; %open figure

for i=1:numPlaces
    
    %display gerberas who had the lowest score on this statistic
    subplot(3,numPlaces,i);
    imshow([imgDir labels{lowestIdxs(i)}   '.png']);
    title([placeLabels{i} ' L: ' num2str(sprintf('%3.4g',lowestValues(i)))]);
    axis image;
    xlabel(labels{lowestIdxs(i)})

    %display gerberas who had mid-level (or close to median) scores on this statistic
    subplot(3,numPlaces,numPlaces+i);
    imshow([imgDir labels{midIdxs(i)} '.png']);
    title(['M: ' num2str(sprintf('%3.4g',midValues(i)))]);
    axis image;
    xlabel(labels{midIdxs(i)})
    
    %display gerberas who had the highest score on this statistic
    subplot(3,numPlaces,3*numPlaces+1-i);
    imshow([imgDir labels{highestIdxs(i)} '.png']);
    title([placeLabels{i} ' H: ' num2str(sprintf('%3.4g',highestValues(i)))]);
    axis image;
    xlabel(labels{highestIdxs(i)})
    
end

%creates title for whole figure
ax=axes('Units','Normal','Position',[.075 .09 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on');
title(statLabel,'Interpreter','none');
%h=get(ax,'Title'); %axis handle

%some tweaking of overall image properties
%set(gcf, 'Position', [1, 41, 1920, 963]); %set size on screen
%set(gcf, 'Position', [1, 41, 1620, 563]); %set size on screen
set(gcf, 'Position', [1, 41, 1400, 570]); %set size on screen
%set(gcf, 'Position', [1, 1, 1200, 520]); %set size on screen
saveas(f,[saveDir statLabel '.png']);
close;

end