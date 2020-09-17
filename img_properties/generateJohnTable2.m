function generateJohnTable2

fileDir = [fileparts(mfilename('fullpath')) filesep];
load([fileDir 'RHS_UPOV_Lab_colours.mat']);

header = {'UPOV','Description','R (RGB)','G (RGB)','B (RGB)','H (HSV)','S (HSV)','V (HSV)','L* (L*a*b*)','a* (L*a*b*)','b* (L*a*b*)'};

sorted_UPOV = sortrows(lookuptable,2);

%separate values and put into matrix (not cell)
UPOV   = cell2mat(sorted_UPOV(:,2));
LAB    = cell2mat(sorted_UPOV(:,3:5));

%get mean of each set of UPOV vales
for i=unique(UPOV)'
    vals = find(UPOV==i);
    meanLAB(i,:) = round(mean(LAB(vals,:))); 
end
    
%get one answer per UPOV colour
UPOV   = unique(UPOV);
meanLAB = meanLAB; %(already done)
labels = unique(strtrim(sorted_UPOV(:,6)),'stable');

rgb = lab2rgb(meanLAB);

rgb_bounded = rgb;
rgb_bounded(rgb_bounded>1)=1; %anything above 1 becomes 1
rgb_bounded(rgb_bounded<0)=0; %anything below 0 becomes 0

hsv = rgb2hsv(rgb_bounded);

rgb_int = uint8(rgb*255);
hsv_int = uint8(hsv*255);

%reorder based on euclidean distance
ord = [1,7,69,8,13,70,71,72,47,50,53,57,11,10,56,30,27,25,19,66,62,16,5,2,3,14,15,9,4,12,58,54,51,48,45,49,52,55,46,43,38,41,40,39,60,65,59,63,22,26,29,31,32,42,37,28,36,35,34,33,64,23,24,21,20,18,17,6,67,68,61,73,44];
UPOV = UPOV(ord);
labels = labels(ord);
rgb_int = rgb_int(ord,:);
hsv_int = hsv_int(ord,:);
meanLAB = meanLAB(ord,:);

%generate table for excel
UPOVtable = cell(length(UPOV)+1,11);
UPOVtable(1,:) = header;
UPOVtable(2:end,1) = num2cell(UPOV);
UPOVtable(2:end,2) = labels;
UPOVtable(2:end,3:5) = num2cell(rgb_int);
UPOVtable(2:end,6:8) = num2cell(hsv_int);
UPOVtable(2:end,9:11) = num2cell(meanLAB);

xlswrite('JohnTable2.xlsx',UPOVtable);

end