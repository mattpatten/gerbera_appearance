function plot_histograms_LAB(imgAbbr, gerb, segmentation, radmean_LAB, kmeans_idxs, xCorrVals, yCorrFit, coeffsLAB, corrR, label)


%create/set output directory
outputDir = get_dir(imgAbbr,'output');
LABdistDir = [outputDir 'LABdist' filesep];
if ~exist(LABdistDir,'dir'), mkdir(LABdistDir); end

%set properties
[xSize, ySize, ~] = size(gerb.RGB);
r = convert2polar([xSize ySize]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Compute some values    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3D histogram (frequency per eccentricity for each L*/a*/b*)
histData = get_3D_HistData(gerb.LAB_masked.whole, r);

%convert from index/eccentricity to circle of that eccentricity
kmeans_circles = zeros(size(r)); %initialize
for s=1:length(kmeans_idxs)
    kmeans_circles = kmeans_circles + double(r==kmeans_idxs(s)); %get circle for any eccentricities identified
end
seg_circles = zeros(size(r)) + double(r==segmentation(1)) + double(r==segmentation(2)); %get circle for any eccentricities identified

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         PLOT! :-)         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
rows = 4; cols = 6;

%display flower segments along top row
subplot(rows,cols,1); imshow(applyMask(gerb.RGB,gerb.mask.whole));
subplot(rows,cols,2); imshow(applyMask(gerb.RGB,gerb.mask.disk));
subplot(rows,cols,3); imshow(applyMask(gerb.RGB,gerb.mask.trans));
subplot(rows,cols,4); imshow(applyMask(gerb.RGB,gerb.mask.ray));
subplot(rows,cols,5); imshow(imoverlay(imoverlay(applyMask(gerb.RGB,gerb.mask.whole),kmeans_circles,'g'),seg_circles,'m')); %ray florets with segmentation

%draw histogram of each segment for LUMINANCE channel
draw_histogram_subplot(rows,cols,7, 'Histogram: L',[0 100],gerb.LAB_masked.whole(:,:,1));
draw_histogram_subplot(rows,cols,8, 'Histogram: L',[0 100],gerb.LAB_masked.disk(:,:,1));
draw_histogram_subplot(rows,cols,9, 'Histogram: L',[0 100],gerb.LAB_masked.trans(:,:,1));
draw_histogram_subplot(rows,cols,10,'Histogram: L',[0 100],gerb.LAB_masked.ray(:,:,1));

%draw histogram of each segment for A channel
draw_histogram_subplot(rows,cols,13,'Histogram: A',[-50 100],gerb.LAB_masked.whole(:,:,2));
draw_histogram_subplot(rows,cols,14,'Histogram: A',[-50 100],gerb.LAB_masked.disk(:,:,2));
draw_histogram_subplot(rows,cols,15,'Histogram: A',[-50 100],gerb.LAB_masked.trans(:,:,2));
draw_histogram_subplot(rows,cols,16,'Histogram: A',[-50 100],gerb.LAB_masked.ray(:,:,2));

%draw histogram of each segment for B channel
draw_histogram_subplot(rows,cols,19,'Histogram: B',[-50 100],gerb.LAB_masked.whole(:,:,3));
draw_histogram_subplot(rows,cols,20,'Histogram: B',[-50 100],gerb.LAB_masked.disk(:,:,3));
draw_histogram_subplot(rows,cols,21,'Histogram: B',[-50 100],gerb.LAB_masked.trans(:,:,3));
draw_histogram_subplot(rows,cols,22,'Histogram: B',[-50 100],gerb.LAB_masked.ray(:,:,3));

%Draw 3D histogram plots
subplot(rows,cols,12); hist3(histData{1},'edges',{0:10:(xSize/2)   0:10:100},'CDataMode','auto','FaceColor','interp'); xlabel('Ecc'); xticks(0:50:150); ylabel('L*'); yticks(0:50:100);   zlabel('Freq');
subplot(rows,cols,18); hist3(histData{2},'edges',{0:10:(xSize/2) -50:15:100},'CDataMode','auto','FaceColor','interp'); xlabel('Ecc'); xticks(0:50:150); ylabel('a*'); yticks(-50:50:100); zlabel('Freq');
subplot(rows,cols,24); hist3(histData{3},'edges',{0:10:(xSize/2) -50:15:100},'CDataMode','auto','FaceColor','interp'); xlabel('Ecc'); xticks(0:50:150); ylabel('b*'); yticks(-50:50:100); zlabel('Freq');


%% Radial plots

%Luminance
yLimVals = [min(radmean_LAB(:,1)) max(radmean_LAB(:,1))];
subplot(rows,cols,11); hold all;
plot(xCorrVals{1}                ,yCorrFit{1,1},'c-','LineWidth',2); %linear fit (disk)
plot(xCorrVals{2}+segmentation(1),yCorrFit{1,2},'c-','LineWidth',2); %linear fit (trans)
plot(xCorrVals{3}+segmentation(2),yCorrFit{1,3},'c-','LineWidth',2); %linear fit (ray)
plot(1:segmentation(3),radmean_LAB(1:segmentation(3),1),'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',3); %plot data
plot([segmentation(1) segmentation(1)],yLimVals,'r'); %vertical line between disk and trans florets
plot([segmentation(2) segmentation(2)],yLimVals,'r'); %vertical line between trans and ray florets
line([kmeans_idxs;    kmeans_idxs],   yLimVals,'Color','g'); %vertical line for intra-segmentation of ray florets
box off; axis tight; xlabel('Eccentricity'); ylabel('Mean L*'); ylim(yLimVals); %axis properties
%title(sprintf('%c',bullseye{1}));
title(sprintf('Ray Fl. R^2 = %.2f,  sl = %.2f',corrR{1,3},coeffsLAB{1,3}(1)));

%A-channel
yLimVals = [min(radmean_LAB(:,2)) max(radmean_LAB(:,2))];
subplot(rows,cols,17); hold all;
plot(xCorrVals{1}                ,yCorrFit{2,1},'c-','LineWidth',2); %linear fit (disk)
plot(xCorrVals{2}+segmentation(1),yCorrFit{2,2},'c-','LineWidth',2); %linear fit (trans)
plot(xCorrVals{3}+segmentation(2),yCorrFit{2,3},'c-','LineWidth',2); %linear fit (ray)
plot(1:segmentation(3),radmean_LAB(1:segmentation(3),2),'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',3); %plot data
plot([segmentation(1) segmentation(1)],yLimVals,'r'); %vertical line between disk and trans florets
plot([segmentation(2) segmentation(2)],yLimVals,'r'); %vertical line between trans and ray florets
line([kmeans_idxs;    kmeans_idxs],   yLimVals,'Color','g'); %vertical line for intra-segmentation of ray florets
box off; axis tight; xlabel('Eccentricity'); ylabel('Mean a*'); ylim(yLimVals); %axis properties
%title(sprintf('%c',bullseye{2}));
title(sprintf('Ray Fl. R^2 = %.2f,  sl = %.2f',corrR{2,3},coeffsLAB{2,3}(1)));

%B-channel
yLimVals = [min(radmean_LAB(:,3)) max(radmean_LAB(:,3))];
subplot(rows,cols,23); hold all;
plot(xCorrVals{1}                ,yCorrFit{3,1},'c-','LineWidth',2); %linear fit (disk)
plot(xCorrVals{2}+segmentation(1),yCorrFit{3,2},'c-','LineWidth',2); %linear fit (trans)
plot(xCorrVals{3}+segmentation(2),yCorrFit{3,3},'c-','LineWidth',2); %linear fit (ray)
plot(1:segmentation(3),radmean_LAB(1:segmentation(3),3),'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor','b','MarkerSize',3); %plot data
plot([segmentation(1) segmentation(1)],yLimVals,'r'); %vertical line between disk and trans florets
plot([segmentation(2) segmentation(2)],yLimVals,'r'); %vertical line between trans and ray florets
line([kmeans_idxs;    kmeans_idxs],   yLimVals,'Color','g'); %vertical line for intra-segmentation of ray florets
box off; axis tight; xlabel('Eccentricity'); ylabel('Mean b*'); ylim(yLimVals); %axis properties
%title(sprintf('%c',bullseye{3}));
title(sprintf('Ray Fl. R^2 = %.2f,  sl = %.2f',corrR{3,3},coeffsLAB{3,3}(1)));

set(gcf, 'Position', [1, 1, 2000, 1700]); %set size on screen
saveas(f,[LABdistDir 'LABdist_' label '.png']);
close;

end



function [histData] = get_3D_HistData(maskedImg, r)

[xSize, ~, nColours] = size(maskedImg);

for cc=1:nColours
    histData{cc} = [];
    tmp = maskedImg(:,:,cc);
    for rad=1:(xSize/2)
        vals = tmp(r==rad);
        histData{cc} = [histData{cc}; repmat(rad,length(vals),1) vals];
    end
    histData{cc} = histData{cc}(all(~isnan(histData{cc}),2),:); %remove any NaNs
end
end


function draw_histogram_subplot(figr,figc,figi,figtitle,figEdges,data)

subplot(figr,figc,figi); hold all; title(figtitle); xlim(figEdges); box off;
histogram(reshape(data, numel(data), 1),'BinWidth',1);
%plot([nanmean(data,'all') nanmean(data,'all')],[0 max(histcounts(data,'BinWidth',1))],'r','LineWidth',2);
plot([nanmean(nanmean(data)) nanmean(nanmean(data))],[0 max(histcounts(data,'BinWidth',1))],'r','LineWidth',2);
xlabel('Intensity'); ylabel('Freq');

end