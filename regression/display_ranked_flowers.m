function display_ranked_flowers(imgAbbr,session)

%attribute = 'both';
transpose = false;
nClusters = 1;
cluster_num = 1;
use_zscore = 0;

%get/create directories
[outputDir, dataDir, imgDir] = get_dir(imgAbbr,'output','data','img');
rankDir = [outputDir 'Rankings' filesep];
if ~exist(rankDir,'dir'), mkdir(rankDir); end

if transpose
    load([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_transpose_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat']); %load aesthetic ratings
else
    load([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat']); %load aesthetic ratings
end

if use_zscore
    refRatings = zscore(refRatings,[],2);
end

%sort from max to min
[~, idx] = sort(mean(refRatings),'descend');
sorted_flowerIDs = flwrs(idx);
sorted_ratings   = refRatings(:,idx);
score = mean(sorted_ratings);

[numSubjects, numFlowers] = size(sorted_ratings);

disp(['Number of subjects: ' num2str(numSubjects)]);

f = figure;
row = 5; col = 14;
for i=1:numFlowers
    h(i) = subplot(row,col,i);
    imshow(sprintf('%s%s%04i.png',imgDir,imgAbbr,sorted_flowerIDs(i)));
    title(sprintf('#%i (%.2f)',i,score(i)));
    axis image;
end

width = 0.13;
for i=1:numFlowers
    pos(i,:) = [-0.02+mod(i-1,col)*0.07 0.855-width/2-(floor((i-1)/col))*(width+0.04) width width];
    set(h(i),'Position',pos(i,:)); %[0.1300 0.1100 0.7750 0.8150]
end

set(gcf, 'Position', [1, 41, 1828, 817]); %set size on screen
if transpose
    saveas(f,[rankDir 'rankings_' session '_transpose_cl' num2str(nClusters) '_' num2str(cluster_num) '.png']);
    %saveas(f,[rankDir 'rankings_' session '_transpose_cl' num2str(nClusters) '_' num2str(cluster_num) '.fig']);
else
    saveas(f,[rankDir 'rankings_' session '_cl' num2str(nClusters) '_' num2str(cluster_num) '.png']);
    %saveas(f,[rankDir 'rankings_' session '_cl' num2str(nClusters) '_' num2str(cluster_num) '.fig']);
end
close;

end