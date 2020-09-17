function compareModelsBetweenExps


imgAbbr = 'expone';
analysis_name = 'v6';
ref_session = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
remove_ref_IO = 1;
new_session = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
remove_new_IO = 1;

%clusters
nClusters = 1; %how many clusters was the data broken into
cluster_num = 1; %which cluster group do we want to run

%shuffling
shuffling = 0; %whether to perform shuffling. Choose either this or boostrapping, not both.
bootstrapping = 1; %whether to perform bootstrapping. Choose either this or shuffling, not both.
%nIterations = 1000; %number of shuffles/bootstraps to perform when performing shuffling or bootstrapping

%transformations
rankTransform = 0; %whether to leave as raw values (0) or convert to ranks (1)
use_zscore = 1; %whether to use a z-score transformation on subjects on subject data


%%%%%%%%%%%%%%%%%%%%%%
%%    Initialize    %%
%%%%%%%%%%%%%%%%%%%%%%

%save dir
[outputDir, imgDir] = get_dir(imgAbbr,'output','img');
comparisonDir = [outputDir 'Experiment_Comparisons' filesep 'ref_IO' num2str(remove_ref_IO) filesep];
if ~exist(comparisonDir,'dir'), mkdir(comparisonDir); end

%% load untested data
flowers_to_remove = get_influential_observers(new_session,remove_new_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, new_session, remove_new_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' new_session '.mat']);

%rename key variables
new_data = refRatings;
new_b    = lm.Coefficients.Estimate;
new_btstrp_b = b;
%new_pred_bstrp_l = pred_bstrp_lower;
%new_pred_bstrp_u = pred_bstrp_upper;

if remove_new_IO
    new_IO_ratings = mean(excl_refRatings,1);
    new_IO_labels  = excl_labels;
end


%% load reference data
flowers_to_remove = get_influential_observers(ref_session,remove_ref_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, ref_session, remove_ref_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' ref_session '.mat']);

%rename key variables
ref_data = refRatings;
ref_b    = lm.Coefficients.Estimate;
ref_btstrp_b = b;
ref_pred_bstrp_l = pred_bstrp_lower;
ref_pred_bstrp_u = pred_bstrp_upper;

clear refRatings; %just to make sure we don't use it accidentally


%% Important data
trained_data      = mean(ref_data,1)';
untested_data     = mean(new_data,1)';
model_predictions = lm.Fitted;


%%%%%%%%%%%%%%%%%%%%%%%%
%%    Scatterplots    %%
%%%%%%%%%%%%%%%%%%%%%%%%

%% Residuals
%{
ceilfix = @(x) (ceil(abs(x)*10).*sign(x))/10;

resall = [trained_data-model_predictions; untested_data-model_predictions];
buffer = 0.1;

fitmin = ceilfix(min(model_predictions)-buffer);
fitmax = ceilfix(max(model_predictions)+buffer);
resmin = ceilfix(min(resall)-buffer);
resmax = ceilfix(max(resall)+buffer);


f = figure;

%intra-residuals
subplot(2,2,1); hold all;
line([fitmin fitmax],[0 0],'LineStyle','--');
scatter(model_predictions,trained_data-model_predictions,'x','LineWidth',1);
set(gca,'XLim',[fitmin fitmax],'YLim',[resmin resmax]);
xlabel('Fitted values','Interpreter','None');
ylabel('Residuals','Interpreter','None');
title([ref_session '_residuals_on_' ref_session '_fit'],'Interpreter','None');
set(gca,'linewidth',1);
box off;

%inter-residuals
subplot(2,2,3); hold all;
line([fitmin fitmax],[0 0],'LineStyle','--');
scatter(model_predictions,untested_data-model_predictions,'x','LineWidth',1);
set(gca,'XLim',[fitmin fitmax],'YLim',[resmin resmax]);
xlabel('Fitted values','Interpreter','None');
ylabel('Residuals','Interpreter','None');
title([new_session '_residuals_on_' ref_session '_fit'],'Interpreter','None');
set(gca,'linewidth',1);
box off;

%intra-histogram
subplot(2,2,2);
histogram(trained_data-model_predictions,'BinWidth',0.05);
set(gca,'XLim',[resmin resmax]);
xlabel('Residuals','Interpreter','None');
ylabel('Freq','Interpreter','None');
title([ref_session '_residuals_on_' ref_session '_fit'],'Interpreter','None');
set(gca,'linewidth',1);
box off;

%inter-histogram
subplot(2,2,4);
histogram(untested_data-model_predictions,'BinWidth',0.05);
set(gca,'XLim',[resmin resmax]);
xlabel('Residuals','Interpreter','None');
ylabel('Freq','Interpreter','None');
title([new_session '_residuals_on_' ref_session '_fit'],'Interpreter','None');
set(gca,'linewidth',1);
box off;

set(gcf, 'Position', [1, 41, 1228, 917]); %set size on screen
saveas(f,[comparisonDir 'residuals_' new_session '_IO' num2str(remove_new_IO) '_residuals_on_' ref_session '_IO' num2str(remove_ref_IO) '_fit.png']);
close;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Predict & compare    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[~, sort_idx] = sort(b(2:end),'descend'); %Sort by beta values
%title_text = sprintf('Regression Coefficients - Cluster %d of %d, %s%s%s%s  r²=%.02f,  F=%.02f,  p=%.02f,  errvar=%.02f',...
%    cluster_num,nClusters,shuffling_text{1+shuffling},bootstrap_text{1+bootstrapping},rank_text{1+rankTransform},zscore_text{1+use_zscore},rg.stats(1),rg.stats(2),rg.stats(3),rg.stats(4));
%saveName = [RSAdir 'regression_' shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} 'cl' num2str(nClusters) '_' num2str(cluster_num)];



%find model prediction for flowers - rating doesn't matter, this only uses predictors i.e., ypred = predict(mdl,Xnew)
%model_predictions = predict(lm,T); %this literally is just getting model predictions for the same flowers - we have this, lm.Fitted !

%R-squared via correlation coefficient
tmp = corr([trained_data model_predictions]);
corr_Rsquared = tmp(1,2).^2;

%Matlab central function
r2 = rsquare(trained_data,model_predictions);

disp(['R2 for existing data - ' ref_session ' - (corr^2 / rsquare / lm / regress):']);
disp([corr_Rsquared r2 lm.Rsquared.Ordinary rg.stats(1) ]);

%% On new experiment
tmp = corr([untested_data model_predictions]);
corr_Rsquared = tmp(1,2).^2;

%Matlab central function
r2 = rsquare(untested_data,model_predictions);

disp(['R2 for untested data - ' new_session ' on ' ref_session ' model - (corr^2 / rsquare):']);
disp([corr_Rsquared r2]);

plotModelFit(untested_data,model_predictions,ref_pred_bstrp_l,ref_pred_bstrp_u,r2,[new_session '_IO' num2str(remove_new_IO) '_data_on_' ref_session '_IO' num2str(remove_ref_IO) '_fit'],imgDir,comparisonDir,labels,0,new_session);
plotModelFit(untested_data,model_predictions,ref_pred_bstrp_l,ref_pred_bstrp_u,r2,[new_session '_IO' num2str(remove_new_IO) '_data_on_' ref_session '_IO' num2str(remove_ref_IO) '_fit'],imgDir,comparisonDir,labels,1,new_session);


end