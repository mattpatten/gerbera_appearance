function [bad_sIDs] = flowerpoll_QA(attribute,outputDir,subjectID_lims,numQs,dispFig)

% This function takes the CSV files downloaded from the flowerpoll online database (an online survey that asked about different flower properties on prolific)
% and saves the data to a .mat file to be read into our processing pipeline
%attribute = 'interest'; %are we extracting 'appeal', 'interest' or 'physical' information
% 
% Created by Matt Patten
% Cleaned up in July 2020


%% Parameters
unresponsive_limit = 50; %>= maximum per cent of answers for a single value, before we no longer consider that this person had enough variety in their responses to consider
max_sds_from_mean  = 2; %number of standard deviations from the mean, where any lower it is considered that they did not have enough variance in their answers to be doing the task correctly, or responded just way too quickly.
numAttributes = 12; %appeal, bullseye etc - everything present in the model, whether it is used in the specific analysis or not
numNonDataCols = 14; %subjectID, consent, order, linkID, gender, ... (non-data columns prior to main data)

physical_attributes = {'bullseye','busyness','complexity','depth','pointiness','symmetry'};


%% File I/O
databaseDir = get_dir([],'database'); %where to access the database files from
saveDir = [outputDir 'QA' filesep];
if ~exist(saveDir,'dir'), mkdir(saveDir); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Part 1: QA data file   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open file
QA.filename = 'quality_check.csv';
fid = fopen([databaseDir QA.filename]);

%get headers
QA.raw_headers = textscan(fid,repmat('%s ',1,numNonDataCols+4*numAttributes),1,'delimiter',',','CollectOutput',1);
QA.raw_headers = QA.raw_headers{1};

%get data
QA.raw = textscan(fid,['%f %s %q %s %s %s %f %s %f %f %s %f %f %f ' repmat('%f %f %f %f ',1,numAttributes)],'delimiter',',');

%close file
fclose(fid);

%Extract attribute variables (i.e., remove demographic variables)
property_idxs = (numNonDataCols+1):length(QA.raw);
QA.headers    = QA.raw_headers(property_idxs);
QA.attributes = cell2mat(QA.raw(:,property_idxs));

%define columns
subj_col         = find(ismember(QA.raw_headers,{'subjectID'})); 
unresponsive_col = ismember(QA.raw_headers,{'unresponsive'});

%cut down subjects
first_ppt_num  = subjectID_lims(1);
last_ppt_num   = min([subjectID_lims(2) max(QA.raw{subj_col})]); %min so we can set last value as "too big" and will just take last value as we continue testing and add more ppts
subj_rows      = and(QA.raw{subj_col}>=first_ppt_num, QA.raw{subj_col}<=last_ppt_num); 
QA.attributes  = QA.attributes(subj_rows,:);

nSubjects = sum(subj_rows);

%cut down missing information from questions we didn't ask (i.e., remove nan's row-by-row)
for ss=1:nSubjects
    cleanvers(ss,:) = QA.attributes(ss,~isnan(QA.attributes(ss,:)));
end
QA.attributes = cleanvers;
%QA.headers = QA.headers(any(~isnan(QA.attributes))); %check for nan's by checking if a whole column of data is empty
%QA.attributes = QA.attributes(:,all(~isnan(QA.attributes))); %check for nan's along whole columns


%cut down and isolate relevant demographic information
QA.subjectIDs   = QA.raw{subj_col}(subj_rows);
QA.unresponsive = QA.raw{unresponsive_col}(subj_rows);

%convert min and max response values into a range  (output: range1 / mean1 / sd1 / range2 / mean2 / sd2)
if numQs==1
    QA.final_headers = {'Range P1','Mean P1','SD P1'};
    QA.final = [QA.attributes(:,2)-QA.attributes(:,1), QA.attributes(:,3), QA.attributes(:,4)];
elseif numQs>1
    QA.final_headers = {'Range P1','Mean P1','SD P1','Range P2','Mean P2','SD P2'};
    QA.final = [QA.attributes(:,2)-QA.attributes(:,1), QA.attributes(:,3), QA.attributes(:,4) ...
                QA.attributes(:,6)-QA.attributes(:,5), QA.attributes(:,7), QA.attributes(:,8)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Part 2: Response times   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get headers
RT.filename = 'resp_times.csv';
fid = fopen([databaseDir RT.filename]);
RT.headers = textscan(fid,repmat('%s ',1,6),1,'delimiter',',','CollectOutput',1);
RT.headers = RT.headers{1};

%get data (the old way for my old matlab version)
%RT.data = readmatrix([fileDir filename_resp],'NumHeaderLines',1,'Delimiter',','); %more modern way not available on my home matlab
RT.data = textscan(fid,repmat('%f ',1,6),'delimiter',',');
RT.data = double(cell2mat(RT.data));
fclose(fid);

%cut down data to relevant rows
RT.data = RT.data(subj_rows,:);

%split data
RT.subjectIDs = RT.data(:,1);
RT.data      = RT.data(:,2:end);
RT.headers   = RT.headers(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Part 3: Proportion of each value   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any([strcmpi(attribute,'appeal') strcmpi(attribute,'interest')])

    PROP.filename = ['props_' attribute '.csv'];

    %get headers
    fid = fopen([databaseDir PROP.filename]);
    PROP.headers = textscan(fid,repmat('%s ',1,12),1,'delimiter',',','CollectOutput',1);
    PROP.headers = PROP.headers{1};

    %get data
    PROP.data = textscan(fid,repmat('%f ',1,12),'delimiter',',');
    PROP.data = double(cell2mat(PROP.data));
    fclose(fid);

    %cut down data to relevant rows
    PROP.data = PROP.data(subj_rows,:);

    %split data
    PROP.subjectIDs = PROP.data(:,1);
    PROP.data      = PROP.data(:,2:end);
    PROP.headers   = PROP.headers(2:end);
    
    %which column of the output we want to use to search for unresponsiveness
    PROP.target_col_header = 'prop5';
    
elseif strcmpi(attribute,'physical')
    
    %extract data from QC file seeing as we don't have a prop file for these variables
    PROP.subjectIDs  = QA.subjectIDs;
    PROP.data       = QA.unresponsive;
    PROP.headers    = 'unresponsive';

    %which column of the output we want to use to search for unresponsiveness
    PROP.target_col_header = 'unresponsive';

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Part 4: Analyse      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check subject IDs match
if any(std([QA.subjectIDs RT.subjectIDs PROP.subjectIDs],[],2))
    error('Subject IDs do not match between csv files');
else
    disp(['All subject IDs match. (n=' num2str(nSubjects) ')'])
end

%combine charts
% Range / Mean / SD / Min resp / Median resp / Mean resp / Max resp / Total time / proportion answered as 0,1,2..,5(unresponsive),...10
COMB.data = [QA.final RT.data PROP.data];
COMB.headers = [QA.final_headers RT.headers PROP.headers];

%compute mean and standard deviation limits
COMB.mean = nanmean(COMB.data,1);
COMB.sd   = nanstd(COMB.data,[],1);
COMB.sd_lower = COMB.mean - COMB.sd * max_sds_from_mean;
COMB.sd_upper = COMB.mean + COMB.sd * max_sds_from_mean;

%compute median and min/max using quartiles
COMB.median     = quantile(COMB.data,0.5);
COMB.dist_lower = quantile(COMB.data,0.25) - 1.5*iqr(COMB.data);
COMB.dist_upper = quantile(COMB.data,0.75) + 1.5*iqr(COMB.data);

%get column for unresponsiveness
unresp_col = find(ismember(COMB.headers, PROP.target_col_header));

%put it all together
plot_cols = [1:(5+3*numQs) unresp_col];

%plot
if dispFig, plot_data(saveDir, attribute, COMB, plot_cols, numQs, 'full'); end

%{
%find all those that lie external to 2 s.d. on *ALL* properties
for i=1:size(COMB.data,2)
    outlier_idxs{i}  = [find(COMB.data(:,i)<COMB.sd_lower(i)); find(COMB.data(:,i)>COMB.sd_upper(i))]; 
    outlier_sIDs{i}   = [subjectIDs(COMB.data(:,i)<COMB.sd_lower(i)); subjectIDs(COMB.data(:,i)>COMB.sd_upper(i))];
end
bad_idxs = unique(sort(vertcat(outlier_idxs{1:6}))); %select which columns we want to use for selection here
bad_sIDs = unique(sort(vertcat(outlier_sIDs{1:6}))); %select which columns we want to use for selection here
%}

% Manual override (for first dataset only):
%COMB.sd_lower(3) = 0.666972;
%COMB.sd_lower(6) = 0.735726;
%COMB.sd_lower(8) = 1.998788;


%% display, find and remove outliers

outlier_idxs = [];
outlier_sIDs = []; 

%sd for param 1
sd1_col = find(ismember(COMB.headers,{'SD P1'})); 
outlier_idxs = [outlier_idxs;          find(COMB.data(:,sd1_col)<COMB.sd_lower(sd1_col))];
outlier_sIDs = [outlier_sIDs; QA.subjectIDs(COMB.data(:,sd1_col)<COMB.sd_lower(sd1_col))];
fprintf('\n%0.1f standard deviations below mean for SD parameter 1: %f',max_sds_from_mean, COMB.sd_lower(sd1_col));

%sd for param 2 (if is there)
if numQs==2
    sd2_col = find(ismember(COMB.headers,{'SD P2'})); 
    outlier_idxs = [outlier_idxs;          find(COMB.data(:,sd2_col)<COMB.sd_lower(sd2_col))];
    outlier_sIDs = [outlier_sIDs; QA.subjectIDs(COMB.data(:,sd2_col)<COMB.sd_lower(sd2_col))];
    fprintf('\n%0.1f standard deviations below mean for SD parameter 2: %f',max_sds_from_mean, COMB.sd_lower(sd2_col)); 
end

%median
median_col = find(ismember(COMB.headers,{'median_resp_time'}));
outlier_idxs = [outlier_idxs;          find(COMB.data(:,median_col)<COMB.sd_lower(median_col))];
outlier_sIDs = [outlier_sIDs; QA.subjectIDs(COMB.data(:,median_col)<COMB.sd_lower(median_col))];
fprintf('\n%0.1f standard deviations below mean for median response time: %f',max_sds_from_mean, COMB.sd_lower(median_col));
fprintf('\nMaximum unresponsive rate: %d\n', unresponsive_limit);

%unresponsive
if ismember(attribute,{'appeal','interest'})

    %check each prop column
    for i=0:10
        prop_col = find(ismember(COMB.headers,{sprintf('prop%i',i)}));
        outlier_idxs = [outlier_idxs;          find(COMB.data(:,prop_col)>=unresponsive_limit)];
        outlier_sIDs = [outlier_sIDs; QA.subjectIDs(COMB.data(:,prop_col)>=unresponsive_limit)];
    end

elseif strcmpi(attribute,'physical')

    %find column for unresponsive as calculated by Postgres
    unresp_col = find(ismember(COMB.headers,{'unresponsive'}));
    outlier_idxs = [outlier_idxs;          find(COMB.data(:,unresp_col)>=unresponsive_limit)];
    outlier_sIDs = [outlier_sIDs; QA.subjectIDs(COMB.data(:,unresp_col)>=unresponsive_limit)];

end

bad_idxs = unique(sort(outlier_idxs));
bad_sIDs = unique(sort(outlier_sIDs));

COMB.data(bad_idxs,:) = []; %delete rows

if dispFig, plot_data(saveDir, attribute, COMB, plot_cols, numQs, 'cut'); end %re-plot after cuts

disp(['Subject IDs that did not fit criteria (' num2str(length(bad_sIDs)) '): ']);
disp(bad_sIDs);

end


function plot_data(saveDir, attribute, COMB, cols, numQs, vv)

%hack :-(
numQs_text = {'','both_'};
if strcmpi(attribute,'physical'), numQs_text = {'',''}; end

descriptive_titles = repmat({'Range','Mean','S.D.'},1,numQs);
titles = {descriptive_titles{:},'Min Resp (s)','Median Resp (s)','Mean Resp (s)','Max Resp (s)','Total Time (min)','Unresponsiveness (%)'};
lower_limit = [repmat([-0.5   0   0],1,numQs) 0   0  0   0  0  0];
upper_limit = [repmat([9.5  10 4.5],1,numQs) 10  15 20 120 25 60];
bin_width   = [repmat([1   0.5 0.2],1,numQs) 0.5 0.5  1   5  1  3];
%ylim_upper  = repmat(200,1,9); %[250 200 150 250 150 150 140 100 150 150 150 150];

f = figure; 
for i=1:numel(cols)
    this_data = COMB.data(:,cols(i));
    this_data(this_data>upper_limit(i))=upper_limit(i)+bin_width(i)/2;
    subplot(2+numQs,3,i); histogram(this_data,'BinLimits',[lower_limit(i) upper_limit(i)+bin_width(i)],'BinWidth',bin_width(i)); box off; title(titles{i}); ylabel('Freq');
    xlim([lower_limit(i) upper_limit(i)+bin_width(i)]);
    %ylim([0 ylim_upper(i)]);
    yl = ylim;
    line([COMB.median(cols(i))      COMB.median(cols(i))],  yl,'Color','b',                   'LineStyle','-', 'LineWidth',0.5); %median
    line([COMB.mean(cols(i))        COMB.mean(cols(i))],    yl,'Color','red',                 'LineStyle','-', 'LineWidth',1);   %mean
    line([COMB.sd_lower(cols(i))    COMB.sd_lower(cols(i))],yl,'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',1.5); %lower sd
    line([COMB.sd_upper(cols(i))    COMB.sd_upper(cols(i))],yl,'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',1.5); %upper sd
    %line([COMB.dist_lower(cols(i)) COMB.dist_lower(cols(i))],yl,'Color','b','LineStyle',':','LineWidth',1);   %lower dist (quartiles)
    %line([COMB.dist_upper(cols(i)) COMB.dist_upper(cols(i))],yl,'Color','b','LineStyle',':','LineWidth',1);   %upper dist (quartiles)
    clear this_data;
end

set(gcf, 'Position', [100, 50, 1600, 900]); %set size on screen
saveas(f,[saveDir 'histograms_' numQs_text{numQs} attribute '_bounded_' vv '.png']);
%saveas(f,[saveDir 'histograms_' attribute '_bounded_' vv '.fig']);
close;
end