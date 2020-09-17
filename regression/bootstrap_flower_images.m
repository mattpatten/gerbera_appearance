function bootstrap_flower_images(imgAbbr,nClusters,cluster_num)


%% Parameters
numBootstraps = 1000;
%nClusters = 3;
%cluster_num = 1;
shuffling = 0;
attribute_name = 'appeal';

rng('shuffle');

%get/create directories
[outputDir, dataDir, imgDir] = get_dir(imgAbbr,'output','data','img');
bootstrapDir = [outputDir 'Bootstraps' filesep];
if ~exist(bootstrapDir,'dir'), mkdir(bootstrapDir); end

if shuffling 
    shuffling_label = 'shuffled_';
    title_text = 'Shuffled distribution of flower rankings';
else
    shuffling_label = '';
    title_text = 'Bootstraps of flower ranking';
end

%load whole dataset and sort from max to min to have a fixed ordering of flowers
load([dataDir 'refRDMdata_' imgAbbr '_' attribute_name '_REFERENCE_SET.mat'],'refRatings'); %load properties data (candidate RDMs)
[~, descending_idxs{1}] = sort(mean(refRatings),'descend');
clear refRatings;

%load properties table to sort via every other attribute
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels'); %load properties data (candidate RDMs)
for i=1:length(header)
    [~, descending_idxs{i+1}] = sort(imgProperties(:,i),'descend');
end
order_header = [{'appeal'} header];

%load dataset we want to work with and bootstrap
load([dataDir 'refRDMdata_' imgAbbr '_' attribute_name '_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat']); %load properties data (candidate RDMs)

[numSubjects, numFlowers] = size(refRatings);

%Bootstrap/shuffle the values for these participants, and get the revelant percentiles on either end
for b=1:numBootstraps
    
    clear idxs;

    %shuffle ratings/rankings
    if shuffling
        for i=1:size(refRatings,1)
            %samp(i,:) = refRatings(i,randperm(size(refRatings, 2))); %ratings
            samp(i,:) = randperm(size(refRatings, 2)); %ranking
        end
    else %bootstrapping
        idxs = randsample(numSubjects,numSubjects,true);
        samp = refRatings(idxs,:);
    end
    [~,temp] = sort(mean(samp),'descend');
    [~, rank(b,:)] = sort(temp);
end

for att=1:length(descending_idxs)

    %Re-sort (to whole dataset, max to min)
    sorted_flowerIDs = flwrs(descending_idxs{att});
    sorted_rank      = rank(:,descending_idxs{att});
    
    percentile_lower  = prctile(sorted_rank, 2.5);
    percentile_median = prctile(sorted_rank, 50);
    percentile_upper  = prctile(sorted_rank, 97.5);
    
    %% Draw figure
    
    f = figure;
    hold all;
    
    errorbar(1:numFlowers,percentile_median,percentile_median-percentile_lower,percentile_upper-percentile_median, ...
        'LineStyle','None','LineWidth',1,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','None','Marker','o','MarkerSize',3,'CapSize',3); %parameters
    
    %Make it look nice
    title(sprintf('%s (n=%d, m=%d) - Sorted by %s - Cluster %d of %d',title_text,numSubjects,numBootstraps,order_header{att},cluster_num,nClusters),'Interpreter','None');
    set(gca,'XLim',[0 numFlowers+1],'XTick',1:numFlowers,'XTickLabels',[]);
    set(gca,'YLim',[0 71],'YTick',[1 10:10:70]); %,'YDir','reverse');
    xtickangle(45);
    %xlabel('Flower ID');
    ylabel('Ranking');
    box off; %axis square;
    
    set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
    flwr_width=0.035;
    for i=1:numFlowers
        axes('pos',[0.13-flwr_width/2+i/71*0.775 0.11-flwr_width-flwr_width*mod(i+1,2) flwr_width flwr_width]) %[0.1300 0.1100 0.7750 0.8150]
        imshow(sprintf('%s%s%04i.png',imgDir,imgAbbr,sorted_flowerIDs(i)));
    end
    
    saveas(f,[bootstrapDir 'bootstraps_' shuffling_label '_cl' num2str(nClusters) '_' num2str(cluster_num) '_n' num2str(numSubjects) '_b' num2str(numBootstraps) '_' order_header{att} '.png']);
    close;
end
end