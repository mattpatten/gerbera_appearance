function generate_regression_table


imgAbbr = 'expone';
analysis_name = 'v6';

session1 = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
remove_IO_1 = 1;

session2 = 'appeal'; %used to load the right dataset appeal/interest/both_appeal/both_interest
remove_IO_2 = 0;
        
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

%get reference model
outputDir = get_dir(imgAbbr,'output');
comparisonDir = [outputDir 'Experiment_Comparisons' filesep]; %for histogram only
if ~exist(comparisonDir,'dir'), mkdir(comparisonDir); end
savename = [session1 '_IF' num2str(remove_IO_1) '_to_' session2 '_IF' num2str(remove_IO_2)];

%% load second dataset first (so it will save in directory of first session)
flowers_to_remove = get_influential_observers(session2,remove_IO_2);
RSAdir = get_RSAdir(outputDir, analysis_name, session2, remove_IO_2, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);

load([RSAdir  'modeldata_' imgAbbr '_' session2 '.mat']);
T2 = Te;


%% load first dataset
flowers_to_remove = get_influential_observers(session1, remove_IO_1);
RSAdir = get_RSAdir(outputDir, analysis_name, session1, remove_IO_1, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num);

load([RSAdir  'modeldata_' imgAbbr '_' session1 '.mat']);
T1 = Te;

clear T Te; %to avoid accidentally using these later


T1.Properties.VariableNames = cellfun(@(x) [x '1'],T1.Properties.VariableNames,'UniformOutput',0); %change table headings (so they're not the same between the tables)
T2.Properties.VariableNames = cellfun(@(x) [x '2'],T2.Properties.VariableNames,'UniformOutput',0); 

%add blank row if necessary
emptyAtt = find(~ismember(T1.Attribute1,T2.Attribute2)); %Are there missing entries?
if ~isempty(emptyAtt)
    T2 = [T2(1:(emptyAtt-1),:); cell2table({T1.Attribute1(emptyAtt),NaN,''},'VariableNames',T2.Properties.VariableNames); T2((emptyAtt):end,:)];
end

%{
%% Add rank change to IO table
[~,pos] = ismember(T2.Attribute2,T1.Attribute1);  %where does the IO attribute fall on the ALL table?
rankChange = compose('%+d', pos - (1:size(T2,1))');   %calculate rank change
rankChange(strcmp('+0',rankChange))={'0'};              %replace +0 with 0
rankChange(emptyAtt) = {'---'};

T2 = [T2 rankChange]; %concatenate tables
T2.Properties.VariableNames{end} = 'Rank_Change'; %add heading to new column
%}

%% Combine tables 
%T_full = [T1 T2];
%writetable(T_full,[RSAdir 'regress_compare_full_' session1 '_IO' num2str(remove_IO_1) '_' session2 '_IO' num2str(remove_IO_2) '.xls']);

[T_mini, ord] = join(T1,T2,'LeftKeys',1,'RightKeys',1);
if ~isempty(emptyAtt)
    ord(ord==emptyAtt) = NaN;
    ord(ord>emptyAtt) = ord(ord>emptyAtt)-1; 
end
T_mini = [T_mini num2cell(ord)];
T_mini.Properties.VariableNames{end} = 'Rank'; %add heading to new column
T_mini = T_mini(:,[1 2 3 4 10 6 7]);


%% Beta change

%difference
beta_change = T_mini.Beta_Value2 - T_mini.Beta_Value1;
sd_limit = std(beta_change) .* 2;
big_changes = abs(beta_change) > sd_limit;
asterisks = cell(nAttributes,1);
asterisks(big_changes) = {'*'};

%add to table
%T_mini = [T_mini num2cell(beta_change) asterisks];
%T_mini.Properties.VariableNames((end-1):end) = {'beta_change','beta_change_lt_alpha'}; %add heading to new column
T_mini = [T_mini asterisks];
T_mini.Properties.VariableNames(end) = {'beta_change_lt_2'}; %add heading to new column

create_histogram(beta_change, comparisonDir, savename);


%% Write to excel
T_mini.Attribute1 = rename_headers(T_mini.Attribute1)'; %make labels more presentable
writetable(T_mini,[RSAdir 'regress_compare_' session1 '_IO' num2str(remove_IO_1) '_' session2 '_IO' num2str(remove_IO_2) '.xls']);


%repeat for perceptual and computed properties
%remove_atts = attributes_with_multicollinearity(session,remove_influential_observers); %% Sub-script listing attributes to be removed %%


end



function create_histogram(beta_change,saveDir,saveName)

f = figure;
hold all;
histogram(beta_change,'BinWidth',0.01);
set(gca,'XLim',[-0.2 0.2]);
set(gca,'YTick',0:1:10);
xlabel('Beta change','Interpreter','None');
ylabel('Freq','Interpreter','None');
title(saveName,'Interpreter','None');
box off;
saveas(f,[saveDir 'hist_' saveName '.png']);
close;

end

%{
function RSAdir = get_RSAdir(outputDir, analysis_name, session, remove_influential_observers, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num)

shuffling_text   = {'','shuffled_'};
bootstrap_text   = {'','btstrp_'};
rank_text        = {'','rank_'};
zscore_text      = {'','zscore_'};
influential_text = {'',['_not' sprintf('%d',flowers_to_remove)]};

RSAdir = [outputDir 'RSA' filesep analysis_name '_' session influential_text{1+remove_influential_observers} '_' ...
          shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} ...
          num2str(nClusters) 'cl_' num2str(cluster_num) filesep];

end
%}
 
