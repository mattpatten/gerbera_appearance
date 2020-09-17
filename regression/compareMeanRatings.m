function compareMeanRatings

imgAbbr = 'expone';
analysis_name = 'v6';
ref_session = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
new_session = 'both_appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest

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

remove_IO = 0; %a constant - displays all flowers and IO flowers on all analyses, so no need to ever do this any other way

%save dir
[outputDir, imgDir] = get_dir(imgAbbr,'output','img');
comparisonDir = [outputDir 'Experiment_Comparisons' filesep];
if ~exist(comparisonDir,'dir'), mkdir(comparisonDir); end

%% load untested data

flowers_to_remove = get_influential_observers(new_session,remove_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, new_session, remove_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' new_session '.mat']);

%rename key variables
new_data = refRatings;

%% load reference data
flowers_to_remove = get_influential_observers(ref_session,remove_IO);
RSAdir = get_RSAdir(outputDir, analysis_name, ref_session, remove_IO, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);
load([RSAdir  'modeldata_' imgAbbr '_' ref_session '.mat']);

%rename key variables
ref_data = refRatings;

clear refRatings; %just to make sure we don't use it accidentally


%%%%%%%%%%%%%%%%%%%%%%%%
%%    Scatterplots    %%
%%%%%%%%%%%%%%%%%%%%%%%%

%% Mean Ratings

%get regular and influential observers from dataset
ref_IOs = get_influential_observers(ref_session,1);
main_ref_flwrs = ~ismember(1:size(ref_data,2),ref_IOs); %cut off necessary positions/ranks
IO_ref_flwrs = ~main_ref_flwrs;

new_IOs = get_influential_observers(new_session,1);
main_new_flwrs = ~ismember(1:size(new_data,2),new_IOs); %cut off necessary positions/ranks
IO_new_flwrs = ~main_new_flwrs;

main_flwrs = and(main_new_flwrs,main_ref_flwrs); %

%compute standard errors
new_stderr = std(new_data)/sqrt(size(new_data,1))';
ref_stderr = std(ref_data)/sqrt(size(ref_data,1))';

f = figure;
zmin = -1.75;
zmax = 1.75;
legend_names = {'Line of best fit','Regular flowers','Influential flowers'};
hold all;

%line of best fit
xsteps = linspace(zmin,zmax,200);
polycoef = polyfit(mean(ref_data)',mean(new_data)',1); 
[coef, ~, ~, ~, rating_stats] = regress(mean(new_data)',[mean(ref_data)' ones(size(ref_data,2),1)]);
fitvals = polyval(coef,xsteps); 
plot(xsteps,fitvals,'Color',rgb('LightGrey'),'LineStyle','-','LineWidth',1);

%scatterplot
%scatter(mean(ref_data)',mean(new_data)','LineWidth',1);
errorbar(mean(ref_data(:,main_flwrs))',mean(new_data(:,main_flwrs))',... %data values
         new_stderr(main_flwrs),new_stderr(main_flwrs),ref_stderr(main_flwrs),ref_stderr(main_flwrs),... %error bars
         'Color',rgb('steelblue'),'LineStyle','None','LineWidth',1,'CapSize',0); %formatting

%influential flowers - ref
errorbar(mean(ref_data(:,IO_ref_flwrs))',mean(new_data(:,IO_ref_flwrs))',... %data values
         new_stderr(IO_ref_flwrs),new_stderr(IO_ref_flwrs),ref_stderr(IO_ref_flwrs),ref_stderr(IO_ref_flwrs),... %error bars
         'Color',rgb('orange'),'LineStyle','None','LineWidth',1,'CapSize',0); %formatting

if any(IO_new_flwrs~=IO_ref_flwrs)
    %influential flowers - new
    errorbar(mean(ref_data(:,IO_new_flwrs))',mean(new_data(:,IO_new_flwrs))',... %data values
             new_stderr(IO_new_flwrs),new_stderr(IO_new_flwrs),ref_stderr(IO_new_flwrs),ref_stderr(IO_new_flwrs),... %error bars
             'Color',rgb('forestgreen'),'LineStyle','None','LineWidth',1,'CapSize',0); %formatting
    legend_names = [legend_names{1:2}, {['Influential observers - ' ref_session ' data']}, {['Influential observers - ' new_session ' data']}];
end
     
legend(legend_names,'Location','SouthEast','Interpreter','None');

%Make it look nice
title_text = sprintf('Mean Ratings:  %s & %s,   r²=%.02f,  F=%.02f,  p=%.02f,  errvar=%.02f',ref_session,new_session,rating_stats(1),rating_stats(2),rating_stats(3),rating_stats(4));
title(title_text,'Interpreter','None');
set(gca,'XLim',[zmin zmax],'YLim',[zmin zmax]);
xlabel([ref_session ' (z-score)'],'Interpreter','None');
ylabel([new_session  ' (z-score)'],'Interpreter','None');
set(gca,'linewidth',1);
box off;
axis square;
set(gcf, 'Position', [1, 41, 800, 717]); %set size on screen
saveas(f,[comparisonDir 'meanRatings_' ref_session '_and_' new_session '.png']);
close;

end