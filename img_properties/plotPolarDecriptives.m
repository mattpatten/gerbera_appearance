function plotPolarDecriptives(img, imgAbbr, label)

% This performs functions and displays plots for the descriptive statistics on each colour channel of the polar-transformed image.
%
% Inputs:
%    img   - A struct (.display and .analysis) for a series of images.
%    imgAbbr - 
%    label - Strings providing a name/label for each flower.
%
% Output:
%     A series of graphs analysing the polar transformed images.
%
% Created by Matt Patten in February, 2019.


% Initialize
outputDir = get_dir(imgAbbr,'output');
polarDir = [outputDir filesep 'polarAnalysis' filesep]; 
if ~exist(polarDir,'dir'), mkdir(polarDir); end
nColours = size(img.display,3); %get image properties

f = figure;
rows = 5; cols = 3;
set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
cc_label = {'L*','a*','b*'};
cc_ylims = {[0 100],  [-50 100],  [-50 100]};
cc_ysteps = {0:10:100, -50:10:100, -50:10:100};
subplot(rows,cols,1); imshow(img.display); title('Original');
subplot(rows,cols,2); imshow(lab2rgb(img.analysis)); title('Polar Image');

for cc=1:nColours
    
    %get relevant data
    flwr = img.analysis(:,:,cc);
    
    %mean
    flwrMean{cc} = mean(flwr,1,'omitnan');
    subplot(rows,cols,3+cc);
    plot(1:size(flwr,2),flwrMean{cc},'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',2);
    box off; axis tight; xlabel('Theta'); ylabel(['Mean ' cc_label{cc}]); ylim(cc_ylims{cc});
    
    %mode
    flwrMode{cc} = mode(round(flwr),1);
    subplot(rows,cols,6+cc);
    plot(1:size(flwr,2),flwrMode{cc},'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',2);
    box off; axis tight; xlabel('Theta'); ylabel(['Mode ' cc_label{cc}]); ylim(cc_ylims{cc});
    
    %std
    flwrStd{cc} = std(flwr,[],1,'omitnan'); %std
    subplot(rows,cols,9+cc);
    plot(1:size(flwr,2),flwrStd{cc},'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',2);
    box off; axis tight; xlabel('Theta'); ylabel(['Std ' cc_label{cc}]); ylim([0 30]);
    
    %3d histogram
    flwrHistData{cc} = [reshape(repmat(1:size(flwr,2),size(flwr,1),1),numel(flwr),1) reshape(flwr,numel(flwr),1)]; %one col for theta, one col for the values (so each val is paired with its theta)
    flwrHistData{cc} = flwrHistData{cc}(all(~isnan(flwrHistData{cc}),2),:); %remove any NaNs
    subplot(rows,cols,12+cc);
    hist3(flwrHistData{cc},'edges',{0:24:size(flwr,2) cc_ysteps{cc}},'CDataMode','auto','FaceColor','interp');
    xlabel('Theta (24)'); xticks(0:100:350); ylabel([cc_label{cc} ' (10)']); yticks(cc_ylims{cc}(1):50:cc_ylims{cc}(2)); zlabel('Freq');
end

set(gcf, 'Position', [1, 1, 1000, 1600]); %set size on screen
saveas(f,[polarDir 'polar_' label '.png']);
close;
    
end