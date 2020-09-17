function read_RHS_colours

%loads up xls file containing RHS and CIE Lab colour values.
%sorts it nicely

%load and sort RHS to L*a*b* reference data
fileDir = [fileparts(mfilename('fullpath')) filesep];
[cielab, txt_values] = xlsread([fileDir 'RHS_UPOV_Lab_colours.xlsx'],'CIE Lab');
RHS_lab = txt_values(3:end,1);

%load and sort RHS to UPOV reference data
[~, ~, allData] = xlsread([fileDir 'RHS_UPOV_Lab_colours.xlsx'],'RHS by UPOV');
UPOV = allData(2:end,1);
RHS_upov = allData(2:end,2);
English = allData(2:end,3);

%get indexes of matches
[~, idx_lab, idx_upov] = intersect(RHS_lab,RHS_upov,'stable');

%cut down arrays to only ones that exist in both arrays
RHS_labels     = RHS_upov(idx_upov);
UPOV_labels    = cell2mat(UPOV(idx_upov));
colour_labels  = English(idx_upov);

RHS_labels_lab = RHS_lab(idx_lab); %not needed, but good to have there just to check it duplicates RHS_labels
lab_values     = cielab(idx_lab,:);

%create overall lookup table
lookuptable = [RHS_labels num2cell(UPOV_labels) num2cell(lab_values) colour_labels];
lookuptable_header = {'RHS','UPOV','CIE L*a*b*','Decription'};

%save data
save([fileDir 'RHS_UPOV_Lab_colours.mat'],'RHS_labels','UPOV_labels','colour_labels','lab_values','lookuptable','lookuptable_header');


%Generate a plot of the colorspace (i.e., a 3D plot of the colours available)

rgb_val = lab2rgb(cielab); %convert to rgb
rgb_val(rgb_val<0)=0; %deal with values outside of boundaries
rgb_val(rgb_val>1)=1; %deal with values outside of boundaries

figure;
hold all; grid on;
for i=1:length(rgb_val) %plot each listed colour
    plot3(rgb_val(i,1),rgb_val(i,2),rgb_val(i,3),'LineStyle','None','Marker','o','MarkerSize',6,'MarkerEdgeColor','None','MarkerFaceColor',rgb_val(i,:));
end
xlabel('R'); ylabel('G'); zlabel('B'); title('RHS endorsed colours');

end