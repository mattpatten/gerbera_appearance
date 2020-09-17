function combine_leaveoneout(imgAbbr,nClusters,cluster_num)

numFlowers = 70;
numBootstraps = 1000;

%get/create directories
[dataDir, outputDir] = get_dir(imgAbbr,'data','output');
loadDir = [outputDir 'RSA_leaveoneout_FWE_' num2str(nClusters) 'cl_' num2str(cluster_num) filesep];
%if ~exist(RSAdir,'dir'), mkdir(RSAdir); end

load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels'); %load properties data (candidate RDMs)
[numFlowers, numAttributes] = size(imgProperties);

mean_data = NaN(numAttributes,numFlowers); %pre-allocate
for i=1:numFlowers
    
    %load data and labels
    load([loadDir 'excluding_' num2str(i) filesep 'RSA_output.mat']); %load properties data (candidate RDMs)
    descending_headers = stats.orderedCandidateRDMnames;

    %get original order
    [~, original_idxs] = ismember(header,descending_headers);

    if ~all(strcmp(header',descending_headers(original_idxs)')) %check re-ordering is correct with original array
        error('Not correct order. Debug required.');
    end
    
    %put array back into original order, get its mean (the blue bars) and add it as a 
    %row (each row is missing a single flower in the analysis)
    mean_data(:,i) = mean(stats.candRelatedness_r(:,original_idxs));
  
end

%bootstrapping
bootstrapped_data = NaN(numAttributes,numBootstraps); %pre-allocate
for b=1:numBootstraps 
    idxs = randsample(numFlowers,numFlowers,true);
    bootstrapped_data(:,b) = mean(mean_data(:,idxs),2);
end

percentile_lower  = prctile(bootstrapped_data, 2.5, 2);
%percentile_median = prctile(bootstrapped_data, 50,  2);
percentile_upper  = prctile(bootstrapped_data, 97.5,2);
final_data = mean(mean_data,2);

%% Draw figure

%Sort by mean
[~, sort_idx] = sort(final_data,'descend');
%this_header = header{sort_idx};

f = figure;
hold all;

bar(1:numAttributes,final_data(sort_idx),'FaceColor',[68/225 131/225 149/225],'EdgeColor','None');
errorbar(1:numAttributes,final_data(sort_idx),final_data(sort_idx)-percentile_lower(sort_idx),percentile_upper(sort_idx)-final_data(sort_idx), ...
    'LineStyle','None','LineWidth',1,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','None','Marker','None','MarkerSize',3,'CapSize',3); %parameters

%Make it look nice
title(sprintf('Leave one flower out analysis (b=%d) - Cluster %d of %d',numBootstraps,cluster_num,nClusters),'Interpreter','None');
set(gca,'XLim',[0 numAttributes+1],'XTick',1:numAttributes,'XTickLabels',header(sort_idx),'TickLabelInterpreter','None');
%set(gca,'YLim',[0 71],'YTick',[1 10:10:70],'YDir','reverse');
xtickangle(45);
%xlabel('Flower ID');
set(gca,'ylim',[-0.05 0.25],'ytick',-0.05:0.05:0.25);       
ylabel('RDM Correlation');
box off; %axis square;

set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[loadDir 'bootstraps_cl' num2str(nClusters) '_' num2str(cluster_num) '.png']);
close;

end