function display_rated_flowers(imgAbbr,session,reference_session)

%% Parameters
%session = 'appeal';
use_zscore = 1;
use_ranktransform = 0;
nClusters = 1;
cluster_num = 1;

%get/create directories
[outputDir, dataDir, imgDir] = get_dir(imgAbbr,'output','data','img');
ratingDir = [outputDir 'Ratings' filesep];
if ~exist(ratingDir,'dir'), mkdir(ratingDir); end


%load whole dataset and sort from max to min to have a fixed ordering of flowers
load([dataDir 'flwrpoll_ratings_' imgAbbr '_' reference_session '_REFERENCE_SET.mat']); %load aesthetic ratings

% Apply z-score if necessary
if use_zscore, refRatings = zscore(refRatings,[],2); end

[~, sorted_idxs] = sort(mean(refRatings),'descend');
clear refRatings;

%load dataset we want to work with and bootstrap
load([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat']); %load aesthetic ratings

[numSubjects, numFlowers] = size(refRatings);

% Apply z-score if necessary
if use_zscore, refRatings = zscore(refRatings,[],2); end

%Re-sort (to whole dataset, max to min)
flwrs = flwrs(sorted_idxs);
refRatings = refRatings(:,sorted_idxs);

mean_rating = mean(refRatings);
se_rating    = std(refRatings)/sqrt(numSubjects);


%% Draw figure

f = figure;
hold all;

errorbar(1:numFlowers,mean_rating,se_rating, ...
    'LineStyle','None','LineWidth',1,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','None','Marker','o','MarkerSize',3,'CapSize',3); %parameters

%Make it look nice
title(sprintf('Sorted flower ratings (n=%d) - Cluster %d of %d',numSubjects,cluster_num,nClusters),'Interpreter','None');
set(gca,'XLim',[0 numFlowers+1],'XTick',1:numFlowers,'XTickLabels',[]);
if use_zscore
    set(gca,'YLim',[-2 1]);
elseif use_ranktransform
    set(gca,'YLim',[1 70],'YDir','reverse');
else
    set(gca,'YLim',[1.8 8.2],'YTick',1:8);
end
ylabel('Rating');
box off; %axis square;

set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
flwr_width=0.035;
for i=1:numFlowers
    axes('pos',[0.13-flwr_width/2+i/71*0.775 0.11-flwr_width-flwr_width*mod(i+1,2) flwr_width flwr_width]) %[0.1300 0.1100 0.7750 0.8150]
    imshow(sprintf('%s%s%04i.png',imgDir,imgAbbr,flwrs(i)));
end

saveas(f,[ratingDir 'rating_REF_' reference_session '_rate_' session  '_cl' num2str(nClusters) '_' num2str(cluster_num) '_n' num2str(numSubjects) '.png']);
close;

%disp('Flower rankings: ');
%disp([(1:numFlowers)' flwrs]);

end