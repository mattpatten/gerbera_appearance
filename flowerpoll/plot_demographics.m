function plot_demographics(loadName,saveName,sIDs)

%open file
fid = fopen(loadName);

%get headers
headers = textscan(fid,repmat('%s ',1,11),1,'delimiter',',','CollectOutput',1);
headers = headers{1};

%get data
demographics = textscan(fid,'%d %s %q %s %s %s %d %s %d %d %s','delimiter',',');

%close file
fclose(fid);

%split data
subjectID     = demographics{1};
gender        = categorical(demographics{6}); gender = renamecats(gender,{'Female','Male','Other','Prefer not to say'});
age           = demographics{7};
country       = categorical(demographics{8});
exp_hort      = demographics{9};
exp_fl_des    = demographics{10};
purchase_freq = categorical(demographics{11},{'Bi-weekly','Weekly','Monthly','Quarterly','Never'}); %re-order when input



source = categorical(demographics{4},{'PROLIFIC','SONA','""'}); 
source = renamecats(source,{'Prolific','Sona','Other'}); 
source(isundefined(source)) = 'Other';


%QA - remove bad participants
ppts_to_remove = subjectID(~ismember(subjectID,sIDs)); %anything removed from QA above
keep_idxs = ~ismember(subjectID,ppts_to_remove);

%QA
subjectID     = subjectID(keep_idxs);
source        = source(keep_idxs);
gender        = gender(keep_idxs);
age           = age(keep_idxs);
country       = country(keep_idxs);
exp_hort      = exp_hort(keep_idxs);
exp_fl_des    = exp_fl_des(keep_idxs);
purchase_freq = purchase_freq(keep_idxs);


clear fid;
%save([dataDir 'flwrpoll_demographics_' imgAbbr '_' attribute '_' num2str(nClusters) 'cl_' num2str(clusterNum) '.mat']);


f = figure;
subplot(2,4,1);   histogram(gender);        box off; ylabel('Freq'); title('Gender'); %xtickangle(45);
subplot(2,4,2);   histogram(age);           box off; ylabel('Freq'); title('Age'); set(gca,'XTick',10:10:70);
subplot(2,4,3:4); histogram(country);       box off; ylabel('Freq'); title('Country'); 
subplot(2,4,5);   histogram(source);        box off; ylabel('Freq'); title('Source'); 
subplot(2,4,6);   histogram(exp_hort);      box off; ylabel('Freq'); title('Expertise in Horticulture'); 
subplot(2,4,7);   histogram(exp_fl_des);    box off; ylabel('Freq'); title('Expertise in Floral Design'); 
subplot(2,4,8);   histogram(purchase_freq); box off; ylabel('Freq'); title('Purchase Frequency');

set(gcf, 'Position', [1, 1, 1400, 800]); %set size on screen
saveas(f,saveName);
close;

end