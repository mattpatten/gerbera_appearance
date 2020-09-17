function plot_correlations2(imgAbbr,nClusters,cluster_num)

%idx=40;
%for i=setxor(idx,36:41), plot_correlations('expone',idx,i); end
%Parameters
%nClusters = 1;   %how many clusters was the data broken into
%cluster_num = 1; %which cluster group do we want to run

binWidth = 0.01; %for histogram

%load label and property details from labels file
[dataDir, outputDir] = get_dir(imgAbbr,'data','output');
saveDir = [outputDir 'correlations_and_covariates_' num2str(nClusters) 'cl_' num2str(cluster_num) filesep];
if ~exist(saveDir,'dir'), mkdir(saveDir); end %if directory doesn't exist, create it

%load image properties table
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels');

%load ratings data
load([dataDir 'refRDMdata_' imgAbbr '_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat']); %load properties data (reference RDMs)

%b = {labels{find(and(imgProperties(:,25)>0.94,imgProperties(:,25)<0.95))}}
%show_image_subset('buy',{labels{find(and(imgProperties(:,12)>20,imgProperties(:,12)<21))}})

%threshold = 0.3;
%titlename=header{attribute_row};

[numFlowers, numAttributes] = size(imgProperties);
[numSubjects, numFlowers]   = size(refRatings);

%convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
%rank_rating = convert_to_ranks(refRatings,'descend');
rank_rating = refRatings; %no rank


subjRDM = NaN(numFlowers,numFlowers,numSubjects);
for subj=1:numSubjects
    subjRDM(:,:,subj) = squareform(pdist(rank_rating(subj,:)')); %convert from vector to square matrix / RDMrefRatings
end



%find correlations
for att=1:numAttributes
    
    %convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
    rank_attribute = imgProperties(:,att); %no rank
    %rank_attribute = convert_to_ranks(imgProperties(:,att)','descend')';
    candRDM = squareform(pdist(rank_attribute)); %convert from vector to square matrix / RDM

    for subj=1:numSubjects
        subj_data = subjRDM(:,:,subj);
        allcorr(subj) = corr(squareform(subj_data)',squareform(candRDM)','type','Spearman','rows','pairwise'); %calc
        %subj_data = rank_rating(subj,:)';
        %allcorr(subj) = corr(subj_data,rank_attribute,'type','Spearman'); %calc
    end

    %plot
    f = figure;
    histogram(allcorr,'BinEdges',-1:binWidth:1);
    title(['subj-to-cand: ' header{att}],'Interpreter','None');
    ylabel('Freq');
    xlim([-1 1]);
    box off;
    saveas(f,[saveDir header{att}],'png');
    close;
end
%}


%find covariates
for att2=1:numAttributes

    rank_attribute2 = convert_to_ranks(imgProperties(:,att2)','descend')';

    covarDir = [saveDir header{att2} filesep];
    if ~exist(covarDir,'dir'), mkdir(covarDir); end %if directory doesn't exist, create it
    
    for att1=1:numAttributes
        
        %convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
        %rank_attribute = imgProperties(:,att); %no rank
        rank_attribute1 = convert_to_ranks(imgProperties(:,att1)','descend')';
        %candRDM = squareform(pdist(rank_attribute')); %convert from vector to square matrix / RDM
        
        for subj=1:numSubjects
            %subj_data = subjRDM(:,:,subj);
            %allcorr(subj) = corr(subj_data(:),candRDM(:),'type','Spearman'); %calc
            subj_data = rank_rating(subj,:)';
            covar(subj,1) = corr(subj_data,rank_attribute1,'type','Spearman'); %calc
            covar(subj,2) = corr(subj_data,rank_attribute2,'type','Spearman'); %calc
        end

        %plot
        f = figure;
        scatter(covar(:,2),covar(:,1));
        title(['Covariates. r²=' num2str(corr(covar(:,1),covar(:,2)))], 'Interpreter','None');
        xlabel(header{att2},'Interpreter','None');
        ylabel(header{att1},'Interpreter','None');
        xlim([-1 1]); xticks(-1:0.2:1);
        ylim([-1 1]); yticks(-1:0.2:1);
        box off;
        saveas(f,[covarDir header{att1}],'png');
        close;
    end
end


end