function plot_correlations(imgAbbr,row1,row2)

%idx=40;
%for i=setxor(idx,36:41), plot_correlations('expone',idx,i); end

%load label and property details from labels file
[dataDir, outputDir] = get_dir(imgAbbr,'data','output');
saveDir = [outputDir 'Correlations' filesep];
if ~exist(saveDir,'dir'), mkdir(saveDir); end %if directory doesn't exist, create it

%load image properties table
load([dataDir 'properties_table_' imgAbbr '.mat'],'header','imgProperties','labels');

%b = {labels{find(and(imgProperties(:,25)>0.94,imgProperties(:,25)<0.95))}}
%show_image_subset('buy',{labels{find(and(imgProperties(:,12)>20,imgProperties(:,12)<21))}})

threshold = 0.3;
%row1=36; 
name1=header{row1};
%row2=38; 
name2=header{row2};
titlename=[header{row1} '_vs_' header{row2}];

%find highly correlated functions
allcorr = corr(imgProperties); %calc

%plot
f1 = figure;
imagesc(allcorr); 
colormap('parula'); 
colorbar; 
axis square; 
title('All correlations');
saveas(f1,[saveDir 'All_correlations.png']);
close;

[r, c] = find(allcorr>threshold); %find entries above threshold
highly_correlated = sortrows([r(r~=c) c(r~=c)]); %remove ones down the middle column
highly_correlated = highly_correlated(highly_correlated(:,1)<highly_correlated(:,2),:); %removes double entries [3 6] and [6 3]

disp(' ');
disp(['Parameters that have a correlation above ' num2str(threshold) ': ']);
for i=1:size(highly_correlated,1)
    disp([header{highly_correlated(i,1)} ' & ' header{highly_correlated(i,2)}]);
end

[r, c] = find(allcorr<-threshold); %find entries above threshold
negatively_correlated = sortrows([r(r~=c) c(r~=c)]); %remove ones down the middle column
negatively_correlated = negatively_correlated(negatively_correlated(:,1)<negatively_correlated(:,2),:); %removes double entries [3 6] and [6 3]

disp(' ');
disp(['Parameters that have a correlation below -' num2str(threshold) ': ']);
for i=1:size(negatively_correlated,1)
    disp([header{negatively_correlated(i,1)} ' & ' header{negatively_correlated(i,2)}]);
end
disp(' ');

%scatterplot between two parameters
f2 = figure; 
scatter(imgProperties(:,row1),imgProperties(:,row2)); 
axis square;
xlabel(name1,'Interpreter','None'); 
ylabel(name2,'Interpreter','None'); 
title([titlename ': ' num2str(corr(imgProperties(:,row1),imgProperties(:,row2)))],'Interpreter','None');
saveas(f2,[saveDir titlename '.png']);
close;
end