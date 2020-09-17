function sort_pref_order

% This function loads the exeriment presentation order as specified in the file provided in the root directory
% and re-orders the labelling for the results of the experiment, such that the most preferred stimulus is
% listed first.
%
% Created by Matt Patten in 2018

%load data
rootDir = get_dir([],'root');
load([rootDir 'flowerPrefDataNormed.mat']); %load preference ratings
load([rootDir 'labels_pres.mat']); %load image names in order they were presented

%rename variable
labels_pres = labels; %change the ambiguous 'labels' variable to make it clear that the flowers were presented in this order for the experiment
clear labels; 

%remove missing data (any with negative values)
pptMissingVals = find(min(flowerData,[],2)<0);
flowerData(pptMissingVals,:) = [];

%sort
[meanPrefData, idx] = sort(mean(flowerData),'descend');
sortedFlowerData = flowerData(:,idx);
labels = labels_pres(idx);

%save
save([rootDir 'labels_pref.mat'],'labels','meanPrefData','sortedFlowerData');

end