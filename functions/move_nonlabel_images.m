function move_nonlabel_images(folderName, imgAbbr)

%get image directory
imgDir  = [fileparts(mfilename('fullpath')) filesep 'stimuli' filesep folderName filesep];
maskDir = [imgDir 'masks' filesep];

remDir     = [imgDir 'rem' filesep];
remMaskDir = [remDir 'masks' filesep];

%load (hopefully) cut-down labels file
load([imgDir 'labels_' imgAbbr '.mat']);

%get filenames from directory
a = dir([imgDir imgAbbr '*.png']);
for i=1:length(a)
    filenames{i} = a(i).name;
end

%remove the .png at the end of the filenames
filenames = cellfun(@(x) extractBefore(x,'.'), filenames, 'UniformOutput', 0);

%find  that match
matches = cell2mat(cellfun(@(x) find(strcmp(x,filenames)),labels,'UniformOutput',0));

%find ones that don't overlap
nonmatches = ~ismember(1:length(filenames),matches);

%get their username
files2move = strcat(filenames(nonmatches), '.png');
masks2move = strcat(filenames(nonmatches),'m.png');

%delete all requested files
cellfun(@(d) movefile([imgDir d],[remDir d],'f'), files2move);
cellfun(@(d) movefile([maskDir d],[remMaskDir d],'f'), masks2move);

end