function read_in_flowerpoll_attributes(imgAbbr)

% This function takes the CSV files downloaded from the flowerpoll online database (an online survey that asked about different flower properties on prolific)
% and saves the data to a .mat file to be read into our processing pipeline
%
% session: Can be session,'appeal','interest','both_appeal', or 'both_interest'
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
use_zscore = 0; 
rankTransform = 0; %converts from rating value per subject to ranking value

%Fairly constant parameters
numProp = 12; %appeal, bullseye, busyness, ...
subjectID_lims = [5 177]; 
attribute_headings = {'bullseye','busyness','complexity','depth','pointiness','symmetry'};
%attribute_headings = {'bullseye','busyness','complexity','depth','petal_quantity','petal_size','petal_variability','pointiness','symmetry','uniqueness'};
numQs = 2;
session = 'physical';
    

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Load all ratings   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Quality Assurance   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove all rows that don't have a rating for our desired characteristics - doesn't do much since we've used subject lims anyway
attribute_cols = ismember(all_ratings.headers,attribute_headings);
all_ratings.clean = all_ratings.raw(any(~isnan(all_ratings.raw(:,attribute_cols)),2),:); %now minimal effect seeing as we've selected participant numbers to encompass these attributes

%determine subjects to remove based on our QA
bad_sIDs = flowerpoll_QA(session, outputDir, subjectID_lims, numQs, 1); %run QA to identify bad participants

%remove bad subjects (based on the QA)
%% Users manually deleted for n=50 attribute collection - DON'T DO DOUBLE THE QA (as s.d. will keep shrinking and deleting)
%all_ratings.clean = all_ratings.clean(~ismember(all_ratings.clean(:,1),bad_sIDs),:); %remove participants

%find any database errors (i.e., situations where flowers presented don't match what was shown to the user)
all_sIDs      = unique(all_ratings.clean(:,1));
all_flowerIDs = unique(all_ratings.clean(:,2));

%check for database errors by comparing values for each subject with the range of flowers in the experiment
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


%%%%%%%%%%%%%%%%%%%%%%
%%   Demographics   %%
%%%%%%%%%%%%%%%%%%%%%%

loadName = [databaseDir 'demographics.csv'];
saveName = [demographicsDir 'demographics_' session '.png'];
plot_demographics(loadName,saveName,all_sIDs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Analyse Rating Data   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjID_col = ismember(all_ratings.headers,'subjectID');
flwrID_col = ismember(all_ratings.headers,'flowerID');
attribute_cols = find(~(subjID_col+flwrID_col));

for att=1:length(attribute_cols)

    %get subject IDs who have responded to this attribute
    subj_with_responses = unique(all_ratings.clean(~isnan(all_ratings.clean(:,attribute_cols(att))),subjID_col));
    
    for ss=1:length(subj_with_responses)
        [~, order_idxs] = sort(all_ratings.clean(all_ratings.clean(:,subjID_col)==subj_with_responses(ss),flwrID_col)); %get flower order for this subject
        ratings_per_subj = all_ratings.clean(all_ratings.clean(:,subjID_col)==subj_with_responses(ss),attribute_cols(att)); %get appeal for this subject
        ratings_per_subj = ratings_per_subj(order_idxs); %put ratings in ascending order (sorted by flower ID)
        poll.data{att}(:,ss) = ratings_per_subj; %add to matrix
    end
end

%extract meaningful details out of dataset
for att=1:length(attribute_cols) %could be converted into cellfun/bsxfun when you can be bothered to do so
    poll.means(att,:) = mean(poll.data{att},2);
    poll.sd(att,:)    = std(poll.data{att},[],2);
    poll.count(att,:) = size(poll.data{att},2);
end
poll.headers = all_ratings.headers(attribute_cols);
poll.flowerIDs = flwrs;

save([dataDir 'flwrpoll_ratings_' imgAbbr '_' session '.mat'],'poll');

end
