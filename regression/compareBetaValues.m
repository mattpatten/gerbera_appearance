function compareBetaValues

imgAbbr = 'expone';
analysis_name = 'v5_sameatt';
ref_session = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
new_session = 'both_interest'; %used to load the right dataset appeal/interest/both_appeal/both_interest

remove_IO = 0;

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
[outputDir] = get_dir(imgAbbr,'output');
comparisonDir = [outputDir 'Experiment_Comparisons' filesep];
if ~exist(comparisonDir,'dir'), mkdir(comparisonDir); end

%% load untested data

flowers_to_remove = get_influential_observers(new_session,remove_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, new_session, remove_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' new_session '.mat']);

%rename key variables
new_b    = lm.Coefficients.Estimate(2:end);
new_btstrp_b = b(:,2:end);
new_header = header;


%% load reference data
flowers_to_remove = get_influential_observers(ref_session,remove_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, ref_session, remove_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' ref_session '.mat']);

%rename key variables
ref_b    = lm.Coefficients.Estimate(2:end);
ref_btstrp_b = b(:,2:end);
ref_header = header;

clear refRatings; %just to make sure we don't use it accidentally


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Gather variables    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nRefBeta = length(ref_b);
nNewBeta = length(new_b);

new_b_prct_lower = new_b - prctile(new_btstrp_b, 2.5)';
new_b_prct_upper = prctile(new_btstrp_b,97.5)' - new_b;
ref_b_prct_lower = ref_b - prctile(ref_btstrp_b, 2.5)';
ref_b_prct_upper = prctile(ref_btstrp_b,97.5)' - ref_b;

if nRefBeta~=nNewBeta
    
    missing_att = find(~ismember(ref_header,new_header)); %find which attribute(s) only exist in the reference
    
    %save it separately; remove from original
    extra_b            = ref_b(missing_att);            ref_b(missing_att)            = [];
    extra_b_prct_lower = ref_b_prct_lower(missing_att); ref_b_prct_lower(missing_att) = [];
    extra_b_prct_upper = ref_b_prct_upper(missing_att); ref_b_prct_upper(missing_att) = [];
    nExtras = length(extra_b);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Generate plots    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
blim = 0.23;
hold all;

%line of best fit
xsteps = linspace(-blim,blim,200);
polycoef = polyfit(ref_b,new_b,1); 
[coef, ~, ~, ~, beta_stats] = regress(new_b,[ref_b ones(size(ref_b,1),1)]);
fitvals = polyval(coef,xsteps); 
plot(xsteps,fitvals,'Color',rgb('Silver'),'LineWidth',1);

%plot scatter
%scatter(ref_b,new_b,'LineWidth',1);
errorbar(ref_b,new_b,new_b_prct_lower,new_b_prct_upper,ref_b_prct_lower,ref_b_prct_upper,...
    'Color',rgb('Blue'),'LineStyle','None','LineWidth',1,'CapSize',0);

%if there's an uneven number of betas - plot the ones that only exist on one attribute along zero for the other
if nRefBeta~=nNewBeta
    errorbar(extra_b,zeros(1,nExtras),zeros(1,nExtras),zeros(1,nExtras),extra_b_prct_lower,extra_b_prct_upper,...
        'Color',rgb('Blue'),'LineStyle','None','LineWidth',1,'CapSize',0);
end

%if comparing appeal and interest, text-label any outlying points
if strcmp(new_session,'both_interest')
    header = rename_headers(header);
    refsteps = polyval(coef,ref_b);
    offcentre = abs(new_b-refsteps)>0.07;
    model_headers = [{'Intercept'}, header]';
    text(ref_b(offcentre)+0.003,new_b(offcentre)+0.008,model_headers(offcentre),'Interpreter','None');
end

%Make it look nice
title_text = sprintf('Beta values:  %s & %s,   r²=%.02f,  F=%.02f,  p=%.02f,  errvar=%.02f',ref_session,new_session,beta_stats(1),beta_stats(2),beta_stats(3),beta_stats(4));
title(title_text,'Interpreter','None');
set(gca,'XLim',[-blim blim],'YLim',[-blim blim]);
xlabel([ref_session ' (beta values)'],'Interpreter','None');
ylabel([new_session  ' (beta values)'],'Interpreter','None');
set(gca,'linewidth',1);
box off;
axis square;
set(gcf, 'Position', [1, 41, 800, 717]); %set size on screen
saveas(f,[comparisonDir 'betavalues_' ref_session '_and_' new_session '_IO' num2str(remove_IO) '.png']);
close;

end