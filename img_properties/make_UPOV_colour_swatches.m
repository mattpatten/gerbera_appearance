function make_UPOV_colour_swatches

%generates colour swatches for UPOV chart

fileDir = [fileparts(mfilename('fullpath')) filesep];
load([fileDir 'RHS_UPOV_Lab_colours.mat']);

sorted_UPOV = sortrows(lookuptable,2);

%separate values and put into matrix (not cell)
UPOV   = cell2mat(sorted_UPOV(:,2));
LAB    = cell2mat(sorted_UPOV(:,3:5));

%get mean of each set of UPOV vales
for i=unique(UPOV)'
    vals = find(UPOV==i);
    meanLAB(i,:) = round(mean(LAB(vals,:))); 
end
    
%get one answer per UPOV colour
UPOV   = unique(UPOV);
meanLAB = meanLAB; %(already done)
labels = unique(strtrim(sorted_UPOV(:,6)),'stable');

%sort colours into new order
%[~,idx] = sortrows(meanLAB(:,2)); %by luminance
%UPOV = 	UPOV(idx);
%meanLAB = meanLAB(idx,:);
%labels = labels(idx);


%%{
%choose white, then choose next closest colour
whiteUPOV = 1; 

currentLAB = meanLAB(whiteUPOV,:);

%save data
newUPOV(1)  = UPOV(whiteUPOV);
newLAB(1,:) = meanLAB((UPOV==whiteUPOV),:);
newLabel{1} = labels(whiteUPOV);

%remove from array
UPOV(whiteUPOV)      = [];
meanLAB(whiteUPOV,:) = [];
labels(whiteUPOV)    = [];

for i=2:(length(UPOV)+1)
    [~,distIdx] = min(pdist2(meanLAB,currentLAB)); %compare all values to current value
    currentLAB = meanLAB(distIdx,:); %save current value
    
    %save this value
    newUPOV(i)  = UPOV(distIdx);
    newLAB(i,:) = meanLAB(distIdx,:);
    newLabel{i} = labels(distIdx);
    
    %remove it from array
    UPOV(distIdx)      = [];
    meanLAB(distIdx,:) = [];
    labels(distIdx)    = [];
end

%replace back with original coordinates
UPOV = newUPOV;
meanLAB = newLAB;
labels = newLabel;
%%}


%create cell array of each colour triplet
for cc=1:length(meanLAB)
    img{cc} = repmat(reshape(lab2rgb(meanLAB(cc,:)),1,1,3),100,100);
    %img{cc} = insertText(img{cc},[1 1],num2str(UPOV(cc)),'FontSize',12,'BoxColor','white','BoxOpacity',0.4,'TextColor','black'); %display upov number
    img{cc} = insertText(img{cc},[1 85],labels{cc},'FontSize',10,'BoxColor','white','BoxOpacity',0.4,'TextColor','black'); %display label
end

%tile images
f = montage(img);
title('UPOV colours');

set(gcf, 'Position', [1, 1, 1080, 1080]); %set size on screen
saveas(f,[fileDir,'UPOV_swatches.png']);

end