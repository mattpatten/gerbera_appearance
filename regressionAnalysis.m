function regressionAnalysis(session,remove_influential_observers)

% Performs a regression analysis based on our table of image properties and 
% our behavioural results from the flowerpoll survey.
%
% Inputs:
%    imgAbbr     - An abbrevation linking to the dataset of flowers to process.
%
% Created by Matt Patten
% Created on 28/5/2020


%% Parameters
%labels
imgAbbr = 'expone';
%session = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
analysis_name = 'v6';
%remove_influential_observers = 0;

%clusters
nClusters = 1; %how many clusters was the data broken into
cluster_num = 1; %which cluster group do we want to run

%shuffling
shuffling = 0; %whether to perform shuffling. Choose either this or boostrapping, not both.
bootstrapping = 1; %whether to perform bootstrapping. Choose either this or shuffling, not both.
nIterations = 1000; %number of shuffles/bootstraps to perform when performing shuffling or bootstrapping

%transformations
rankTransform = 0; %whether to leave as raw values (0) or convert to ranks (1)
use_zscore = 1; %whether to use a z-score transformation on subjects on subject data

%lasso regression
lambda_threshold = 0; %for lasso regression and removal of certain paramters (1-100, 0=off)

%correlations
correlation_threshold = 0.5; %maximum level of correlation to allow (identifies them, but doesn't remove them - that's VIF)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Initialize      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

flowers_to_remove = get_influential_observers(session,remove_influential_observers);

shuffling_text = {'','shuffled_'};
bootstrap_text = {'','btstrp_'};
rank_text = {'','rank_'};
zscore_text = {'','zscore_'};
influential_text = {'',['_not' sprintf('%d',flowers_to_remove)]};

if ~or(shuffling,bootstrapping), nIterations = 0; end %only run through once if we're not doing shuffling

%get/create main directory
[outputDir, dataDir, imgDir] = get_dir(imgAbbr,'output','data','img');
RSAdir = [outputDir 'RSA' filesep analysis_name '_' session influential_text{1+remove_influential_observers} '_' ...
          shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} ...
          num2str(nClusters) 'cl_' num2str(cluster_num) filesep];

if ~exist(RSAdir,'dir'), mkdir(RSAdir); end

%load data and labels
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels'); %load properties data (candidate RDMs)
load([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_' num2str(nClusters) 'cl_' num2str(cluster_num) '.mat'],'refRatings'); %load aesthetic ratings

header_disp = rename_headers(header);

%% Convert to z-scores
if use_zscore
    refRatings    = zscore(refRatings,[],2); 
    imgProperties = zscore(imgProperties,[],1);
end


%% Remove influential observations
if remove_influential_observers
    cases_to_include = ~ismember(1:size(refRatings,2),flowers_to_remove); %cut off necessary positions/ranks
    excl_flowers = ~cases_to_include;

    %save excluded cases for later
    excl_labels = labels(excl_flowers);
    excl_refRatings = refRatings(:,excl_flowers);
    excl_imgProperties = imgProperties(excl_flowers,:);
    
    %update
    labels = labels(cases_to_include);
    refRatings = refRatings(:,cases_to_include);
    imgProperties = imgProperties(cases_to_include,:);
end

%correlate one attribute with all others
%n=18; a = corr(imgProperties(:,n),imgProperties); [~, sidx] = sort(abs(a),'ascend'); sidx = sidx(1:(end-1)); disp([num2cell(a(sidx)') header(sidx)']);


%% Remove attributes related to multicollinearity
all_imgProperties = imgProperties; %save to combine with post-removal
all_header = header; 

remove_atts = attributes_with_multicollinearity(session,remove_influential_observers); %% Sub-script listing attributes to be removed %%
imgProperties(:,remove_atts) = [];
header(remove_atts) = []; header_disp(remove_atts) = [];
if remove_influential_observers, excl_imgProperties(:,remove_atts) = []; end


%% Plot correlations heatmaps
not_grad_idxs = ~ismember(1:length(all_header),[23 24]); %this isn't done here because if we split it up it will change the index values for everything else after it
plotCorrHeatmap(corr(all_imgProperties(:,not_grad_idxs)), all_header(not_grad_idxs), [RSAdir 'corr_all_attributes' influential_text{1+remove_influential_observers}]); %all variables
plotCorrHeatmap(corr(imgProperties), header, [RSAdir 'corr_used_attributes' influential_text{1+remove_influential_observers}]); %after removal

corr_lower = tril(corr(all_imgProperties),-1); %get lower triangle
corr_upper = triu(corr(all_imgProperties),0); %get upper triangle, but not main diagonal

%create mask that locates where we need to insert nans
nAttributes = size(all_imgProperties,2);
mask_rem = zeros(nAttributes);
mask_rem(remove_atts,:) = 1; 
mask_rem(:,remove_atts) = 1;
mask = and(triu(ones(nAttributes),0),mask_rem); %only use columns/rows when in the upper triangular

%set specified values in the mask to NaN
corr_upper(mask) = NaN;

%remove gradient attributes
corr_upper = corr_upper(:,not_grad_idxs); corr_upper = corr_upper(not_grad_idxs,:);
corr_lower = corr_lower(:,not_grad_idxs); corr_lower = corr_lower(not_grad_idxs,:); 
all_header = all_header(not_grad_idxs); 

plotCorrHeatmap(corr_upper + corr_lower, all_header, [RSAdir 'corr_cut_upper' influential_text{1+remove_influential_observers}]); %combined version

%{
corr_upper = triu(corr(all_imgProperties),0); %get upper triangle
corr_lower = tril(corr(all_imgProperties),-1); %get lower triangle, but not main diagonal
corr_upper(remove_atts,:) = 0;
corr_upper(:,remove_atts) = 0;
plotCorrHeatmap(corr_upper(:,not_grad_idxs) + corr_lower(:,not_grad_idxs), all_header(not_grad_idxs), [RSAdir 'corr_cut_upper' influential_text{1+remove_influential_observers}]); %combined version
%}

%% Get properties
[nFlowers, nAttributes] = size(imgProperties);
[nSubjects, nFlowers] = size(refRatings);


%% Variance inflation factor (VIF)
vic = vif(imgProperties);

[val,idx] = sort(vic,'ascend');
header_sort = header_disp(idx);

disp('Variance inflation coefficients (sorted, ascending) - threshold of 5 or 10 is recommended: ');
disp([num2cell(val') header_sort']);

vic = vic';

%get index for each attribute
%disp(' ');
%disp([num2cell((1:length(header_disp))') num2cell(vic) header_disp']);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Lasso Regression   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

[lass, lassInfo] = lasso(imgProperties,mean(refRatings,1)','Alpha',1); %lasso regularisation (alpha~=0 is ridge regularisation, inbetween is elastic net)
[ridg, ridgInfo] = lasso(imgProperties,mean(refRatings,1)','Alpha',0.001); %ridge regularisation (alpha~=0 is ridge regularisation, inbetween is elastic net)

nLambdas = size(lass,2);
cmap = jet(nAttributes);

f = figure;
subplot(1,2,1);
hold all;
for att=1:nAttributes
    plot(1:nLambdas,lass(att,:),'Color',cmap(att,:),'LineWidth',1.5);
end
legend(header_disp,'Location','NorthEastOutside','Interpreter','None');
title('Lasso regression','Interpreter','None');
set(gca,'XLim',[1 nLambdas]);
xlabel('Lambda (standardised)','Interpreter','None');
ylabel('Regression Coefficients','Interpreter','None');
set(gca,'linewidth',1);
box off;
axis square;

subplot(1,2,2);
hold all;
for att=1:nAttributes
    plot(1:nLambdas,ridg(att,:),'Color',cmap(att,:),'LineWidth',1.5);
end
legend(header,'Location','NorthEastOutside','Interpreter','None');
title('Ridge regression','Interpreter','None');
set(gca,'XLim',[1 nLambdas]);
xlabel('Lambda (standardised)','Interpreter','None');
ylabel('Regression Coefficients','Interpreter','None');
set(gca,'linewidth',1);
axis square;
box off;

set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[RSAdir 'LassoRidge.png']);
close;

if lambda_threshold>0
    keep_idxs = find(abs(lass(:,lambda_threshold))>0.001);
    disp('Using parameters:');
    disp(header_disp(keep_idxs)');
    
    %cut down parameters
    imgProperties = imgProperties(:,keep_idxs);
    header = header(keep_idxs);
    header_disp = header_disp(keep_idxs);
    [nFlowers, nAttributes] = size(imgProperties);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   High Correlations   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find correlations between attributes
all_corrs = corr(imgProperties);

%find entries that have a correlation above threshold
[r, c] = find(abs(all_corrs)>correlation_threshold); 

%remove ones down the middle column
highly_correlated = sortrows([r(r~=c) c(r~=c)]); 
dif_r = r(r~=c); 
dif_c = c(r~=c);
r = dif_r; clear dif_r;
c = dif_c; clear dif_c;

%removes double entries [3 6] and [6 3]
vals_to_keep = highly_correlated(:,1)<highly_correlated(:,2);
r = r(vals_to_keep);
c = c(vals_to_keep);
highly_correlated = highly_correlated(vals_to_keep,:); 

%sort to display highest correlations first
[~,desc_idx] = sort(abs(all_corrs(sub2ind([nAttributes nAttributes],r,c))),'descend');
highly_correlated = highly_correlated(desc_idx,:);
r = r(desc_idx);
c = c(desc_idx);

disp(' ');
disp(['Highly correlated parameters (>' num2str(correlation_threshold) '): ']);
for i=1:size(highly_correlated,1)
    disp(sprintf('(%.02f) %s & %s',abs(all_corrs(r(i),c(i))),header_disp{highly_correlated(i,1)},header_disp{highly_correlated(i,2)}));
    table_results(i,:) = {abs(all_corrs(r(i),c(i))), header_disp{highly_correlated(i,1)}, header_disp{highly_correlated(i,2)}};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Run bootstrapping / shuffling   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate
b = NaN(nIterations,nAttributes+1); %+ constant/intercept term
fitted_data_bstrps = NaN(nIterations,nFlowers);

for ss=1:nIterations %is set to one if not performing shuffling

    sample_ratings = refRatings; %reset at the start of each loop

    if shuffling
        %shuffle the data and update user
        sample_ratings = Shuffle(sample_ratings')'; 
        disp(['Performing shuffle #' num2str(ss) '...']);

    elseif bootstrapping
        %Bootstrap the values for these participants, and get the revelant percentiles on either end
        bootstrap_idxs = randsample(nSubjects,nSubjects,true); %select values from 1:n, and do this n times, with replacement
        sample_ratings = sample_ratings(bootstrap_idxs,:);
        disp(['Performing bootstrap #' num2str(ss) '...']);
    end %transpose is needed as it does it columnwise

    %% Perform multiple linear regression 
    T = array2table([imgProperties mean(sample_ratings,1)'],'VariableNames',[header, 'appeal_mean'],'RowNames',labels);
    lm = fitlm(T,'linear');
    b(ss,:) = lm.Coefficients.Estimate;
    [rg.b,rg.bint,rg.r,rg.rint,rg.stats] = regress(mean(sample_ratings,1)',[ones(nFlowers,1) imgProperties]); %alternative version
    fitted_data_bstrps(ss,:) = lm.Fitted;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Multiple Linear Regression   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_ratings = mean(refRatings,1)'; %simplify

T = array2table([imgProperties mean_ratings],'VariableNames',[header, session],'RowNames',labels);
lm = fitlm(T,'linear')
[rg.b,rg.bint,rg.r,rg.rint,rg.stats] = regress(mean_ratings,[ones(nFlowers,1) imgProperties]); %alternative version

if bootstrapping
    %get percentiles from bootstraps
    percentile_lower = prctile(b, 2.5);
    mean_val = lm.Coefficients.Estimate';
    percentile_upper = prctile(b,97.5);
    
elseif shuffling

    %sort_vals = lm.Coefficients.Estimate';
    %[~, sort_idx] = sort(sort_vals(2:end),'descend'); %Sort by beta values
    
    %get percentiles from bootstraps
    percentile_lower = prctile(b, 2.5);
    mean_val         = prctile(b,  50);
    percentile_upper = prctile(b,97.5);

    %title_text = sprintf('Regression Coefficients - Cluster %d of %d, %s%s%s%s ',...
    %             cluster_num,nClusters,shuffling_text{1+shuffling},bootstrap_text{1+bootstrapping},rank_text{1+rankTransform},zscore_text{1+use_zscore},rg.stats(1));
    %saveName = [RSAdir 'regression_' shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} 'cl' num2str(nClusters) '_' num2str(cluster_num)];
    %plotRegression(mean_val(2:end), mean_val(2:end)-percentile_lower(2:end), percentile_upper(2:end)-mean_val(2:end), sort_idx, header, title_text, saveName)
    
else
    mean_val = lm.Coefficients.Estimate';
    %[~, sort_idx] = sort(mean_val(2:end),'descend'); %Sort by beta values
    %title_text = sprintf('Regression Coefficients - Cluster %d of %d, %s%s%s%s  r²=%.02f,  F=%.02f,  p=%.02f,  errvar=%.02f',...
    %             cluster_num,nClusters,shuffling_text{1+shuffling},bootstrap_text{1+bootstrapping},rank_text{1+rankTransform},zscore_text{1+use_zscore},rg.stats(1),rg.stats(2),rg.stats(3),rg.stats(4));
    %saveName = [RSAdir 'regression_' shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} 'cl' num2str(nClusters) '_' num2str(cluster_num)];

end

%separate constant/intercept
intercept        = mean_val(1);
percentile_lower = percentile_lower(2:end);
mean_val         = mean_val(2:end);
percentile_upper = percentile_upper(2:end);


[~, sort_idx] = sort(mean_val,'descend'); %Sort by beta values
[~, magsort_idx] = sort(abs(mean_val),'descend'); %Sort by absolute values for beta

% find goodness of fit as you add each regressor
for i=1:length(magsort_idx)

    %[stp.b,stp.bint,stp.r,stp.rint,stp.stats] = regress(mean_ratings,[ones(nFlowers,1) ]); %alternative version
    Tstep = array2table([imgProperties(:,magsort_idx(1:i)) mean_ratings],'VariableNames',[header(magsort_idx(1:i)), session],'RowNames',labels);
    lm_step = fitlm(Tstep,'linear');

    step_obs = eval(['lm_step.Variables.' lm_step.VariableNames{end}]);
    step_pred = lm_step.Fitted;
    
    GoF(i) = sum(((step_obs - step_pred).^2) ./ step_pred); %this is wrong
    Rstep(i) = lm_step.Rsquared.Ordinary;
    Rstep_adj(i) = lm_step.Rsquared.Adjusted;
    %gofFit = goodnessOfFit(step_obs,step_pred,'MSE');
    %line of best fit
    %xsteps = linspace(zmin,zmax,200);
    %polycoef = polyfit(mean(ref_data)',mean(new_data)',1); 
    %[coef, ~, ~, ~, rating_stats] = regress(mean(new_data)',[mean(ref_data)' ones(size(ref_data,2),1)]);
    %fitvals = polyval(coef,xsteps); 

end

Rstep_cell = num2cell(string(cat(2,                            ... %casting
         num2str(Rstep','%.2f'),                                ... %Rstep
         repmat('  (',nAttributes,1),                          ... %space + open bracket sign
         num2str(round(Rstep'/lm.Rsquared.Ordinary.*100),'%d'), ... %percentage of Rstep
         repmat('%)',nAttributes,1))));                            %close bracket

Rstep_adj_cell = num2cell(string((num2str(Rstep_adj','%.2f'))));

     
title_text = sprintf('Regression Coefficients - Cluster %d of %d, %s%s%s%s  r²=%.02f,  F=%.02f,  p=%.02f,  errvar=%.02f',...
    cluster_num,nClusters,shuffling_text{1+shuffling},bootstrap_text{1+bootstrapping},rank_text{1+rankTransform},zscore_text{1+use_zscore},rg.stats(1),rg.stats(2),rg.stats(3),rg.stats(4));
saveName = [RSAdir 'regression_' shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} 'cl' num2str(nClusters) '_' num2str(cluster_num)];

plotRegression(mean_val, mean_val-percentile_lower, percentile_upper-mean_val, sort_idx, header_disp, title_text, saveName)
    

%% Generate Excel table
bCI = compose('[%.2f, %.2f]',percentile_lower',percentile_upper');
Te = cell2table([header_disp(magsort_idx)' num2cell(mean_val(magsort_idx))' bCI(magsort_idx) Rstep_cell Rstep_adj_cell],'VariableNames',{'Attribute','Beta_Value','Bootstrapped_CI','Cumulative_Rsq','Adjusted_Rsq'});
writetable(Te,[saveName '.xls']);


%% Compute for perceptual / objective attributes only

poll_attributes = {'poll_bullseye','poll_busyness','poll_complexity','poll_depth','poll_pointiness','poll_symmetry'};
poll_cols = ismember(header,poll_attributes);

for j=1:2
    perc_imgProperties = imgProperties;
    perc_header = header;
    
    %measure perceptual properties
    perc_imgProperties = perc_imgProperties(:,poll_cols);
    perc_header = perc_header(poll_cols);

    T_perc = array2table([perc_imgProperties mean_ratings],'VariableNames',[perc_header, session],'RowNames',labels);
    lm_perc = fitlm(T_perc,'linear');
    Rperc(j) = lm_perc.Rsquared.Ordinary;
    Rperc_adj(j) = lm_perc.Rsquared.Adjusted;
    nPerceptual(j) = sum(poll_cols);
    poll_cols = ~poll_cols;
end

T_perc = cell2table([{'Perceptual';'Objective'} num2cell(Rperc') num2cell(Rperc_adj') num2cell(nPerceptual')],'VariableNames',{'Analysis','R_Squared_Ord','R_Squared_Adj','Num_Attributes'});
writetable(T_perc,[RSAdir 'perceptual_regression.xls']);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Generate RDMs    %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%generate reference RDMs from *individual* preference psychological data
%this is better because it looks at individual differences between flower preferences, which is then averaged across everybody
refRDMvect = NaN((nFlowers^2-nFlowers)/2,nSubjects);
refRDM = NaN(nFlowers,nFlowers,nSubjects);
for i=1:nSubjects
    refRDMvect(:,i) = pdist(refRatings(i,:)');
    refRDM(:,:,i)   = squareform(refRDMvect(:,i));
end

%generate candidate RDMs from the flower image properties
candRDMvect = NaN((nFlowers^2-nFlowers)/2,nAttributes);
candRDM     = NaN(nFlowers,nFlowers,nAttributes);
for att=1:nAttributes
    candRDMvect(:,att) = pdist(imgProperties(:,att));
    candRDM(:,:,att)   = squareform(candRDMvect(:,att));
end

%normalize - not sure this works like this - do check
%candRDM = candRDM - min(min(candRDM)); %make minimum 0 (likely is already the case)
%candRDM = candRDM / max(max(candRDM)); %make maximum 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   RDM Correlation matrix   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
corrmat = NaN(nSubjects,nSubjects);   %preallocate

for subj=1:nSubjects
    corrmat(subj,:) = corr(refRDMvect(:,subj),refRDMvect,'type','Pearson','rows','pairwise');
end

%for display purposes, remove perfect correlation down middle diagonal
corrmat_disp = corrmat;
corrmat_disp(corrmat_disp>0.999)=NaN;

f = figure;
colormap(colormapFromRSA);
imagesc(corrmat_disp);
title('Subject correlation matrix','Interpreter','None');
set(gca,'XTick',0:100:nSubjects,'YTick',0:100:nSubjects);
xlabel('Subject'); ylabel('Subject');
axis square;
colorbar;
set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[RSAdir 'corrmat.png']);
close;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Display all RDMs   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
f = figure; rows = 7; cols = 6;
all_headers = [session header_disp];
all_RDMs    = cat(3,mean(refRDM,3),candRDM);
for att=1:(nAttributes+1)
    subplot(rows,cols,att);
    colormap(colormapFromRSA);
    imagesc(rankTransform_equalsStayEqual(all_RDMs(:,:,att)));
    title(all_headers{att},'Interpreter','None');
    axis square off;
end
set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[RSAdir 'allRDMs.png']);
close;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Rank-Ratings and Rank-Attribute correlations   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
RDMrankDir = [RSAdir 'RDMrankcor' filesep];
if ~exist(RDMrankDir,'dir'), mkdir(RDMrankDir); end

%meanRatings_rank = convert_to_ranks(mean(refRDMvect,2)','descend')'; %convert RDM values to ranks, averaged across people
meanRatings_rank = convert_to_ranks(mean(refRatings,1),'descend');

for att=1:nAttributes
    
    %meanAtt_rank = convert_to_ranks(candRDMvect(:,att)','descend')';
    meanAtt_rank = convert_to_ranks(imgProperties(:,att)','descend');
    f = figure;
    hold all;
    
    scatter(meanRatings_rank,meanAtt_rank);
    
    %Make it look nice
    set(gca,'Xlim',[1 length(meanRatings_rank)],'YLim',[1 length(meanAtt_rank)]);
    title(['Rank correlation: ' header_disp{att}], 'Interpreter','None');
    xlabel('Subject Rank','Interpreter','None');
    ylabel([header_disp{att} ' rank'],'Interpreter','None');
    axis square;
    box off;
    
    set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
    saveas(f,[RDMrankDir header{att} '.png']);
    close;
end
%}


%%%%%%%%%%%%%%%%%%%%%%%%
%%   Rank transform   %%
%%%%%%%%%%%%%%%%%%%%%%%%

%convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
if rankTransform
    refRatings = convert_to_ranks(refRatings,'descend');
    if ss==1
        for att=1:nAttributes
            %convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
            imgProperties(:,att) = convert_to_ranks(imgProperties(:,att)','descend')';
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Influential Observations   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfluenceDir = [RSAdir 'Influence' filesep];
if ~exist(InfluenceDir,'dir'), mkdir(InfluenceDir); end

dffits_max = 2;
[~,flwr_rank] = sort(lm.Variables{:,end},'descend');
sorted_dffits = lm.Diagnostics.Dffits(flwr_rank);

f = figure; hold all;

%draw horizontal bounds
line([0 nFlowers+1],[ dffits_max  dffits_max],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',2);
line([0 nFlowers+1],[-dffits_max -dffits_max],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',2);

%draw data
plot(1:nFlowers,sorted_dffits,'LineStyle','None','Marker','x','MarkerSize',12,'LineWidth',2,'Color',rgb('red'));
%scatter(1:nFlowers,sorted_dffits,'Marker','x','LineWidth',2,'MarkerEdgeColor',rgb('red'));

%draw vertical line connecting x-axis to outliers
ylim_val = get(gca,'Ylim');
val = find(abs(sorted_dffits)>dffits_max)';
line([val; val],[repmat(ylim_val(1),1,length(val)); sorted_dffits(val)'],'Color',rgb('LightGrey'),'LineWidth',1.5);

%decorate
title(['DFFITS to identify influential observers - ' session],'Interpreter','None');
ylabel('Scaled change in fit');
set(gca,'XLim',[0 nFlowers+1],'XTick',1:nFlowers,'XTickLabels',[]);
set(gca,'linewidth',1);

box off; %axis square;
set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
flwr_width=0.035; %nFlowers/200
for flwr=1:nFlowers
    axes('pos',[0.13-flwr_width/2+flwr/(nFlowers+1)*0.775 0.11-flwr_width-flwr_width*mod(flwr+1,2) flwr_width flwr_width]) %[0.1300 0.1100 0.7750 0.8150]
    imshow(sprintf('%s%s.png',imgDir,labels{flwr_rank(flwr)}));
end

saveas(f,[RSAdir 'DFFITS.png']);
close;

% Matlab standard figures

%influence_plots = {'contour','cookd','covratio','dfbetas','dffits','leverage','s2_i'};
influence_plots = {'dffits','cookd'};
for i=1:length(influence_plots)
    plotDiagnostics(lm,influence_plots{i});
    if strcmpi(influence_plots{i},'dffits')
        %delete Matlab threshold
        children = get(gca, 'children'); 
        delete(children(1));
    end
    f = gcf; box off;
    saveas(f,[InfluenceDir influence_plots{i} '.png']);
    close;
end

%[val,idx] = sort(lm.Diagnostics.CooksDistance,'descend');
[~,idx] = sort(abs(lm.Diagnostics.Dffits),'descend');
disp('DFFITS values (descending): ');
disp([num2cell(lm.Diagnostics.Dffits(idx)) labels(idx)']);


%%%%%%%%%%%%%%%%%%%%%%%
%%     Residuals     %%
%%%%%%%%%%%%%%%%%%%%%%%

ResidualDir = [RSAdir 'Residuals' filesep];
if ~exist(ResidualDir,'dir'), mkdir(ResidualDir); end


%residual_plots = {'caseorder','fitted','histogram','lagged','probability','symmetry'};
residual_plots = {'fitted','histogram'};
for i=1:length(residual_plots)
    plotResiduals(lm,residual_plots{i});
    f = gcf; box off;
    saveas(f,[ResidualDir residual_plots{i} '.png']);
    close;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Simple linear regressions   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
RegressionDir = [RSAdir 'Regression' filesep];
if ~exist(RegressionDir,'dir'), mkdir(RegressionDir); end

for att=1:nAttributes
    
    %[b,bint,r,rint,stats] = regress(mean(refRatings,1)',[ones(nFlowers,1) imgProperties(:,att)]);
    xsteps = linspace(min(imgProperties(:,att)),max(imgProperties(:,att)),100);
    
    f = figure;
    hold all;
    
    scatter(imgProperties(:,att),mean(refRatings,1));
    plot(xsteps,b(1) + b(1+att)*xsteps,'r','LineWidth',2);
    
    %Make it look nice
    title(sprintf('Multiple linear regression:  %s  b=%0.2f',header_disp{att},b(1+att)),'Interpreter','None');
    set(gca,'XLim',[xsteps(1) xsteps(end)]);
    xlabel(header_disp{att},'Interpreter','None');
    ylabel('Appeal (z-score)','Interpreter','None');
    box off;
    axis square;
    set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
    saveas(f,[RegressionDir header{att} '.png']);
    close;
        
    %% Segmented regression
    
    if and(att,~shuffling)
        SegRegDir = [RSAdir 'Segmented Regression' filesep];
        if ~exist(SegRegDir,'dir'), mkdir(SegRegDir); end
    end
    
    slm = slmengine(imgProperties(:,att),mean(refRatings,1)','degree',1,'plot','on','knots',3,'interiorknots','free');
    f = gcf;
    title(header_disp{att},'Interpreter','None');
    xlabel([header_disp{att} ' (z-score)'],'Interpreter','None');
    ylabel('Mean Appeal (z-score)',   'Interpreter','None');
    box off;
    saveas(f,[SegRegDir header{att} '.png']);
    close;

end
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Model evaluation   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

observed_data = eval(['lm.Variables.' lm.VariableNames{end}]);
%predicted_data = lm.Fitted;
pred_bstrp_med   = prctile(fitted_data_bstrps,  50)';
pred_bstrp_lower = pred_bstrp_med - prctile(fitted_data_bstrps, 2.5)';
pred_bstrp_upper = prctile(fitted_data_bstrps,97.5)' - pred_bstrp_med;
plotModelFit(observed_data,pred_bstrp_med,pred_bstrp_lower,pred_bstrp_upper,lm.Rsquared.Ordinary,session,imgDir,RSAdir,labels,0,session);
plotModelFit(observed_data,pred_bstrp_med,pred_bstrp_lower,pred_bstrp_upper,lm.Rsquared.Ordinary,session,imgDir,RSAdir,labels,1,session);


%% Test model predictions from influential observers
%{
if remove_influential_observers
    excludedT = array2table([excl_imgProperties mean(excl_refRatings,1)'],'VariableNames',[header, 'excluded_appeal_mean'],'RowNames',excl_labels);
    excluded_observed_data = excludedT{:,end};
    excluded_predicted_data = predict(lm,excludedT);
    excl_rsquare = rsquare(excluded_observed_data,excluded_predicted_data);
    plotModelFit(excluded_observed_data,excluded_predicted_data,0,0,excl_rsquare,['excl_' session],imgDir,RSAdir,excl_labels,0);
    plotModelFit(excluded_observed_data,excluded_predicted_data,0,0,excl_rsquare,['excl_' session],imgDir,RSAdir,excl_labels,1);
end
%}

%% Examine values for individual attributes between influential and non-influential flowers
%{
if remove_influential_observers

    %set directory and create folder, if necessary
    AttributesDir = [RSAdir 'Attributes' filesep];
    if ~exist(AttributesDir,'dir'), mkdir(AttributesDir); end
    
    attribute_ROI = 'GLCM_Tamura_coarseness';
    attribute_ROI_disp = rename_headers({attribute_ROI});

    idx = ismember(header,attribute_ROI); %find index for this attribute

    %% Single plot
    f = figure;
    hold all;

    xlim_val = [-4 4]; %get(gca,'Xlim');
    ylim_val = [-2 2]; %get(gca,'Xlim');
    xsteps = linspace(xlim_val(1),xlim_val(2),200);

    %lines of best fit
    %influential
    polycoef = polyfit(excl_imgProperties(:,idx),mean(excl_refRatings)',1); 
    [~, ~, ~, ~, rating_stats] = regress(excl_imgProperties(:,idx),[mean(excl_refRatings)' ones(size(excl_refRatings,2),1)]);
    fitvals = polyval(polycoef,xsteps); 
    plot(xsteps,fitvals,'Color',rgb('LightGrey'),'LineStyle','-','LineWidth',1);

    %non-influential
    polycoef = polyfit(imgProperties(:,idx)',mean_ratings',1); 
    [~, ~, ~, ~, rating_stats] = regress(imgProperties(:,idx),[mean_ratings ones(size(mean_ratings,1),1)]);
    fitvals = polyval(polycoef,xsteps); 
    plot(xsteps,fitvals,'Color',rgb('LightGrey'),'LineStyle','-','LineWidth',1);

    %data
    plot(imgProperties(:,idx),mean_ratings,'LineStyle','None','Marker','x','MarkerSize',12,'LineWidth',2,'Color',rgb('red')); %non-influential
    plot(excl_imgProperties(:,idx),mean(excl_refRatings),'LineStyle','None','Marker','x','MarkerSize',12,'LineWidth',2,'Color',rgb('blue')); %influential
    
    %decorate
    set(gca,'Xlim',xlim_val,'YLim',ylim_val);
    title(sprintf('Both sets of flowers r=%.2f',rating_stats(1)),'Interpreter','None');
    xlabel(attribute_ROI_disp,'Interpreter','None');
    ylabel(session,'Interpreter','None');
    set(gca,'linewidth',1);
    axis square;
    box off;
    
    %set(gcf, 'Position', [1, 41, 1228, 517]); %set size on screen
    saveas(f,[AttributesDir attribute_ROI '2.png']);
    close;
        
    %{
    %% Subplots
    f = figure;
    subplot(1,2,1); hold all;
    plot(excl_imgProperties(:,idx),mean(excl_refRatings),'LineStyle','None','Marker','x','MarkerSize',12,'LineWidth',2,'Color',rgb('red'));
    
    %line of best fit
    xlim_val = [-4 4]; %get(gca,'Xlim');
    ylim_val = [-2 2]; %get(gca,'Xlim');
    xsteps = linspace(xlim_val(1),xlim_val(2),200);
    polycoef = polyfit(excl_imgProperties(:,idx),mean(excl_refRatings)',1); 
    [~, ~, ~, ~, rating_stats] = regress(excl_imgProperties(:,idx),[mean(excl_refRatings)' ones(size(excl_refRatings,2),1)]);
    fitvals = polyval(polycoef,xsteps); 
    plot(xsteps,fitvals,'Color',rgb('LightGrey'),'LineStyle','-','LineWidth',1);

    %decorate
    set(gca,'Xlim',xlim_val,'YLim',ylim_val);
    title(sprintf('Influential flowers r=%.2f',rating_stats(1)),'Interpreter','None');
    xlabel(attribute_ROI_disp,'Interpreter','None');
    ylabel(session,'Interpreter','None');
    set(gca,'linewidth',1);
    axis square;
    box off;

    subplot(1,2,2); hold all;
    plot(imgProperties(:,idx),mean_ratings,'LineStyle','None','Marker','x','MarkerSize',12,'LineWidth',2,'Color',rgb('red'));

    %line of best fit
    %xlim_val = get(gca,'Xlim');
    xsteps = linspace(xlim_val(1),xlim_val(2),200);
    polycoef = polyfit(imgProperties(:,idx)',mean_ratings',1); 
    [~, ~, ~, ~, rating_stats] = regress(imgProperties(:,idx),[mean_ratings ones(size(mean_ratings,1),1)]);
    fitvals = polyval(polycoef,xsteps); 
    plot(xsteps,fitvals,'Color',rgb('LightGrey'),'LineStyle','-','LineWidth',1);

    %decorate
    set(gca,'Xlim',xlim_val,'YLim',ylim_val);
    title(sprintf('Non-influential flowers r=%.2f',rating_stats(1)),'Interpreter','None');
    xlabel(attribute_ROI_disp,'Interpreter','None');
    ylabel(session,'Interpreter','None');
    set(gca,'linewidth',1);
    axis square;
    box off;
    
    set(gcf, 'Position', [1, 41, 1228, 517]); %set size on screen
    saveas(f,[AttributesDir attribute_ROI '.png']);
    close;
    %}
end
%}

save([RSAdir  'modeldata_' imgAbbr '_' session '.mat']);

end


%%%%%%%%%%%%%%%%%%%%%%
%%   Run (my) RSA   %%
%%%%%%%%%%%%%%%%%%%%%%

%{
%convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
if rankTransform, refRatings = convert_to_ranks(refRatings,'descend'); end

relatedness = NaN(nSubjects,nAttributes); %intialize
for att=1:nAttributes
    
    %convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
    if rankTransform, imgProperties(:,att) = convert_to_ranks(imgProperties(:,att)','descend')'; end
    
    %The whole RSA - Compare each subject's RDM to the current candidate RDM
    relatedness(:,att) = corr(refRDMvect,candRDMvect(:,att),'type','Spearman','rows','pairwise')';
    
end

%calculate means
mean_relatedness(ss,:)   = mean(relatedness);
stderr_relatedness(ss,:) = std(relatedness)/sqrt(nSubjects);

%Wilcoxon signed rank test
p_values(att) = signrank(mean(refRDMvect,2),candRDMvect(:,att));

end

%Bonferroni FWE testing
asterisks = p_values<(0.05/nAttributes);
%text(test, y_sorted(test)+es(test)+0.01, ['\bf\fontsize{16}',deunderscore(pStringCell{test})], 'Rotation', 0, 'Color', [0 0 0],'HorizontalAlignment','center');

%Sort by mean
[~, sort_idx] = sort(mean(mean_relatedness,1),'descend');
%this_header = header{sort_idx};

save([RSAdir 'MyRSA_output.mat']);
    

title_text = sprintf('(My) RSA - Cluster %d of %d, %s%s%s',cluster_num,nClusters,shuffling_text{1+shuffling},rank_text{1+rankTransform},zscore_text{1+use_zscore});
saveName = [RSAdir 'relatedness_' shuffling_text{1+shuffling} rank_text{1+rankTransform} zscore_text{1+use_zscore} 'cl' num2str(nClusters) '_' num2str(cluster_num)];
if shuffling
    errorbar_lower  = prctile(mean_relatedness,50) - prctile(mean_relatedness, 2.5);
    errorbar_higher = prctile(mean_relatedness, 97.5) - prctile(mean_relatedness,50);
    plotRegression(prctile(mean_relatedness,50), errorbar_lower, errorbar_higher, sort_idx, header, title_text, saveName); 
else
    errorbar_lower = stderr_relatedness;
    errorbar_higher = stderr_relatedness;
    plotRegression(mean(mean_relatedness,1), errorbar_lower, errorbar_higher, sort_idx, header, title_text, saveName); 
end

end
%}

function rankMat=rankTransform_equalsStayEqual(mat)

% transforms the matrix mat by replacing each element by its rank in the
% distribution of all its elements. if scale01 is set the ranked elements
% would be scaled between 0 to 1. One property of this function, as
% indicated by its name, is that equal values would stay equal after
% ranking and scaling. This feature is not present in the other function of
% the toolbox (rankTransform.m). The way the function does this is by
% assigning the mean scaled rank to all the equal entries of the input
% matrix (mat). NaNs are ignored in this process.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

nonNan_LOG=~isnan(mat);
set=mat(nonNan_LOG);

[sortedSet, sortedIs]=sort(set);

rankMat=nan(size(mat));
nonNan_IND=find(nonNan_LOG);
rankMat(nonNan_IND(sortedIs))=1:numel(set);

% scale into [0,1]
rankMat=(rankMat-1)/(numel(set)-1);

%% fix artefactual differences
uniqueValues=unique(mat(nonNan_LOG));

for uniqueValueI=1:numel(uniqueValues)
    cValueEntries_LOG=mat==uniqueValues(uniqueValueI);
    rankMat(cValueEntries_LOG)=mean(rankMat(cValueEntries_LOG));
end

end%function
