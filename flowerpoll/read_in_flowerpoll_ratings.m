function read_in_flowerpoll_ratings(imgAbbr,session)

% This function takes the CSV files downloaded from the flowerpoll online database (an online survey that asked about different flower properties on prolific)
% and saves the data to a .mat file to be read into our processing pipeline
%
% session: Can be 'physical','appeal','interest','both_appeal', or 'both_interest'
% 
% Created by Matt Patten
% Cleaned up in July 2020


% Set up file I/O
[databaseDir, dataDir, outputDir] = get_dir(imgAbbr,'database','data','output');
demographicsDir = [outputDir 'Demographics' filesep];
if ~exist(demographicsDir,'dir'), mkdir(demographicsDir); end
if ~exist(dataDir,'dir'), mkdir(dataDir); end

%dbDir = 'E:\Matt\Postgres\pg_output\';

%Parameters
nClusters = 1;
transpose = false; %regular - cluster participants, transpose - cluster flowers
use_zscore = 0; 
use_corrmat = 0; %find difference between each user and every other user
rankTransform = 0; %converts from rating value per subject to ranking value

%Fairly constant parameters
numProp = 12; %appeal, bullseye, busyness, ...
att_param = 1;
    
if strcmpi(session,'appeal')
    subjectID_lims = [180 1279]; 
    attribute_headings = {'appeal'}; 
    numQs = 1;    
    
elseif strcmpi(session,'interest')
    subjectID_lims = [1285 1374]; 
    attribute_headings = {'interest'}; 
    numQs = 1;
    
elseif strcmpi(session,'both_appeal')
    subjectID_lims = [1376 10000]; 
    attribute_headings = {'appeal','interest'};
    numQs = 2;
    
elseif strcmpi(session,'both_interest')
    subjectID_lims = [1376 10000];
    attribute_headings = {'appeal','interest'}; 
    att_param = 2;
    numQs = 2;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Load all ratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'all_ratings.csv';

%get headers
fid = fopen([databaseDir filename]);
all_ratings.raw_headers = textscan(fid,repmat('%s ',1,4+numProp),1,'delimiter',',','CollectOutput',1);
all_ratings.raw_headers = all_ratings.raw_headers{1};
all_ratings.raw_headers{2} = all_ratings.raw_headers{2}(1:(end-3)); %change subjectID_id to subjectID

%get data
%all_ratings.all_ratings = readmatrix([dbDir filename],'NumHeaderLines',1,'Delimiter',',');

%get data (the old way for my old matlab version)
all_ratings.raw = textscan(fid,['%f %f %f %q ' repmat('%f ',1,numProp)],'delimiter',',');

%get column indexes for parameters we're interested in
relevant_cols = find(ismember(all_ratings.raw_headers,['subjectID' 'flowerID' attribute_headings]));

%update parameters to only include these columns
all_ratings.headers = all_ratings.raw_headers(relevant_cols);
all_ratings.raw     = all_ratings.raw(relevant_cols);
all_ratings.raw     = cell2mat(all_ratings.raw);

%get column location for subjects
subj_col = find(ismember(all_ratings.headers,{'subjectID'})); 

%condense to only the right subjects for this part of the experiment
first_ppt_num = subjectID_lims(1);
last_ppt_num  = min([subjectID_lims(2) max(all_ratings.raw(:,subj_col))]); %min so we can set last value as "too big" and will just take last value as we continue testing and add more ppts
subj_rows = and(all_ratings.raw(:,subj_col)>=first_ppt_num, all_ratings.raw(:,subj_col)<=last_ppt_num); 
all_ratings.raw = all_ratings.raw(subj_rows,:);

fclose(fid);

%if ~strcmpi(attribute,'physical')
%    appeal_col = find(ismember(all_ratings.headers,{'appeal'}),1);
%    interest_col = find(ismember(all_ratings.headers,{'interest'}),1);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: QA on data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove all rows that don't have a rating for our desired characteristic
ratings_col = ismember(all_ratings.headers,attribute_headings);
all_ratings.clean = all_ratings.raw(all(~isnan(all_ratings.raw(:,ratings_col)),2),:); %now minimal effect seeing as we've selected participant numbers to encompass these attributes

%determine subjects to remove based on our QA
bad_sIDs = [];
for i=1:numQs
    cur_bad_sIDs = flowerpoll_QA(attribute_headings{i}, outputDir, subjectID_lims, numQs, 1); %run QA to identify bad participants
    bad_sIDs = [bad_sIDs; cur_bad_sIDs];
end
bad_sIDs = unique(bad_sIDs);

%remove bad subjects (based on the QA)
all_ratings.clean = all_ratings.clean(~ismember(all_ratings.clean(:,1),bad_sIDs),:); %remove participants

%find any database errors (i.e., situations where flowers presented don't match what was shown to the user)
all_sIDs      = unique(all_ratings.clean(:,1));
all_flowerIDs = unique(all_ratings.clean(:,2));

match = [];
for i=1:length(all_sIDs) %for each subject
    flwrs = sort(all_ratings.clean(all_ratings.clean(:,1)==all_sIDs(i),2)); %get flowers presented to this user
    match(i) = isequal(flwrs,all_flowerIDs); %check if same as 70 that are meant to be presented.
end

%Remove any participants where a database writing error occurred
all_ratings.clean = all_ratings.clean(~ismember(all_ratings.clean(:,1),all_sIDs(~match)),:);

%display to user
disp('Database writing errors:');
disp(all_sIDs(~match));

%update what we consider all users now
all_sIDs = unique(all_ratings.clean(:,1)); 

%Sanity check - just double checking the result above
if (sum(all_ratings.clean(:,2)) - sum(repmat(all_flowerIDs',1,numel(all_sIDs))))==0
    disp('HOORAY! After the QAs, the sum of the flower column minus the sum of each set of flowers multiplied by the number of users is zero!');
    disp('This means each user was shown 70 unique flowers and finished the experiment in full.');
else
    disp('ERROR!!! Not everyone has been shown all 70 flowers. Debug necessary.')
end

disp(['Total number of subjects after QA: ' num2str(length(all_sIDs))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Save appeal data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-allocate
refRatings = nan(length(all_sIDs),length(all_flowerIDs)); 

%get right columns
flwrID_col = ismember(all_ratings.headers,'flowerID');
rating_col = ismember(all_ratings.headers,attribute_headings{att_param});

for i=1:length(all_sIDs)

    [~, order_idxs] = sort(all_ratings.clean(all_ratings.clean(:,1)==all_sIDs(i),flwrID_col)); %get flower order for this subject
    ratings_per_subj = all_ratings.clean(all_ratings.clean(:,1)==all_sIDs(i),rating_col); %get rating for this subject
    ratings_per_subj = ratings_per_subj(order_idxs); %put ratings in ascending order (sorted by flower ID)
    refRatings(i,:) = ratings_per_subj; %add to matrix

end

% Apply z-score if necessary
if use_zscore, refRatings = zscore(refRatings,[],2); end

% Convert to ranks
if rankTransform, refRatings = convert_to_ranks(refRatings,'descend'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Correlation matrix   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
[nSubjects, nFlowers] = size(refRatings);
refRDMvect = NaN((nFlowers^2-nFlowers)/2,nSubjects); %pre-allocate

for subj=1:nSubjects
    refRDMvect(:,subj) = pdist(refRatings(subj,:)');
end

for subj=1:nSubjects
    corrmat(subj,:) = corr(refRDMvect(:,subj),refRDMvect,'type','Pearson','rows','pairwise');
end

%for display purposes, remove perfect correlation down middle diagonal
corrmat_disp = corrmat;
corrmat_disp(corrmat_disp>0.999)=NaN;

f = figure;
imagesc(corrmat_disp);
title('Subject correlation matrix','Interpreter','None');
set(gca,'XTick',0:100:nSubjects,'YTick',0:100:nSubjects);
xlabel('Subject'); ylabel('Subject');
axis square;
colorbar;
set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[outputDir 'CorrMat' filesep session '_corrmat.png']);
close;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   K-means clustering & save   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K-means grouping
all_flowers = flwrs;
all_results = refRatings;

if transpose, refRatings = refRatings'; end %transpose to cluster by flowers, not participants
if use_corrmat
    [classes, ~, centroids, cent_dist] = perform_kmeans(corrmat, nClusters, 0);
else
    [classes, ~, centroids, cent_dist] = perform_kmeans(refRatings, nClusters, 0);
end

for i=1:nClusters

    refRatings = refRatings(classes==i,:); %cut down to just this class
    if transpose
        flwrs = flwrs(classes==i);
        refRatings = refRatings';
        save([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_transpose_' num2str(nClusters) 'cl_' num2str(i) '.mat'],'refRatings','flwrs','centroids','cent_dist','all_sIDs');
    else
        save([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '_' num2str(nClusters) 'cl_' num2str(i) '.mat'],'refRatings','flwrs','centroids','cent_dist','all_sIDs');
    end
    refRatings = all_results; %reset back to full results
    flwrs = all_flowers;
    
    %% Plot demographics
    if ~transpose
        loadName = [databaseDir 'demographics.csv'];
        saveName = [demographicsDir 'demographics_' session '_' num2str(nClusters) 'cl_' num2str(i) '.png'];
        cluster_sIDs = all_sIDs(classes==i);
        plot_demographics(loadName,saveName,cluster_sIDs)
    end
end

end