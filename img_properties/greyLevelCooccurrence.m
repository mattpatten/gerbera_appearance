function [stats] = greyLevelCooccurrence(img, p, lbl)

% This measures how often a pixel with the intensity (grey-level) value i occurs in a specific spatial relationship to a pixel with the value j. 
% It is set up to be a horizontal spatial relationship (i.e., a pixel to the right), but can be adjusted to be in any direction and any number of pixels away.
% 
% Algorithm parameters:
%     p.bins         - The number of grey levels to consider. e.g., 32 bins means every 8 grey values are lumped into a single bin (256/32=8 values per bin).
%     p.symmetry     - True/False. Whether we consider an inverse relationship. i.e., [0,1] as well as [1,0]
%     p.offsetValues - Distance to measure between pixels. 1 indicates neighbouring pixels. 4 means there are three pixels in between the two pixels of interest.
%     p.grayLimits   - Image instensity values to scale between (i.e., min and max for greyscale values to consider, the bins will be divided within these values). Default is all: [0 255]
%
% Inputs:
%     img.analysis  - A series of masked, greyscale images in the order of (x, y, 1(=colour), image_index)
%     imgs.display   - For display purposes only, the original colour image.
%     lbls           - What the directory and filenames should be included to describe the specific analysis (.images, .colour, .analysis, .dir)
%
% Created by Matt Patten
% Created in September, 2018


%% File I/O
%[imgDir, maskDir, outputDir] = get_dir(imgAbbr,'img','mask','output');

%labels
stats_Labels = {'Contrast','Correlation','Energy','Homogeneity'}; %Matlab in-built computations for GLCM
stats_YLims = [0 200; -0.5 1; 0 0.35; 0 1]; %y-axis limits for the different statistics.

%pre-allocate
glcm = zeros(p.bins,p.bins,length(p.offsetValues));
stats = zeros(4,length(p.offsetValues));


%% Compute GLCM
warning off; %Matlab shoots back warning that it doesn't consider NaN entries in the GLCM..... many many times
for idx = 1:length(p.offsetValues)
    offset = [0 p.offsetValues(idx)]; %spatial relationship to the pixel being measured: [row col] e.g., [0 1] would be to the right from the reference pixel; [1 1] diagonally up-right
    [glcm(:,:,idx), SI] = graycomatrix(img.analysis,'NumLevels',p.bins,'Offset',offset,'Symmetric',p.symmetry,'GrayLimits',p.grayLimits);
    
    %get statistics and put them in a single array so we don't have to deal with the complicated mess of structs
    current_stats = graycoprops(glcm(:,:,idx));
    stats(1,idx) = current_stats.Contrast;
    stats(2,idx) = current_stats.Correlation;
    stats(3,idx) = current_stats.Energy;
    stats(4,idx) = current_stats.Homogeneity;
end
warning on; %not to silence all warnings - just this one


%% Display results - GLCMs
f = figure;
set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things

%original image
subplot(3,4,1);
imshow(img.display);
title('Original');
axis square; axis off;

%image of greyed flower scaled into n bins and used in matrix
subplot(3,4,2);
imagesc(SI);
title(['Scaled: ' num2str(p.bins) ' bins']);
colormap('gray');
axis square; axis off;

%grey level cooccurrence matrices
graphsToShow = [1 2 4 8 16 32];
for i=1:length(graphsToShow)
    f = subplot(3,4,2+i);
    imagesc(glcm(:,:,graphsToShow(i)));
    title(['Offset: ' num2str(graphsToShow(i))]);
    axis square;
end

%statistics of GLCM (across offsets)
for j=1:size(stats,1)
    f = subplot(3,4,2+length(graphsToShow)+j);
    plot(1:length(p.offsetValues),stats(j,:),'LineStyle','-','LineWidth',1,'Color','b');
    title(stats_Labels{j});
    set(gca,'YLim',stats_YLims(j,:),'XLim',[1 length(p.offsetValues)],'XTick',1:10:length(p.offsetValues)); %'XTickLabel',glcm_offsetValues);
    xlabel('Offset');
    axis square;
end

%some tweaking of overall image properties
set(gcf, 'Position', [1, 1, 1400, 1250]); %set size on screen
saveas(f,[lbl.dir lbl.image '_' lbl.analysis '.png']);
close;

%{
%% Display results - Between-flower statistics
offsetsToShow = [1 2 8];
colours = [0 0 0; 0.5 0.5 0.5; 0.7 0.7 0.7];
lineWidths = [1 0.7 0.5];
markerStyle = {'-','--',':'};
for off=1:length(offsetsToShow)
    legendName{off} = [num2str(offsetsToShow(off)) ' pixels'];
end

f = figure;
for j=1:size(stats,1) %for each statistic
    subplot(2,2,j);
    hold all;
    for off=1:length(offsetsToShow) %put each shown offset on new line
        plot(1:nImages,squeeze(stats(j,offsetsToShow(off),:)),'LineStyle',markerStyle{off},'LineWidth',lineWidths(off),'Color',colours(off,:));
    end
    if j==1 %print legend in first graph only
        legend(legendName,'Location','NorthEast');
    end
    title(stats_Labels{j});
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages); %'XTickLabel',glcm_offsetValues);
    xlabel('Flower ID');
end

%some tweaking of overall image properties
set(gcf, 'Position', [1, 1, 900, 700]); %set size on screen
saveas(f,[lbl.dir 'stats_' lbl.analysis '.png']);
saveas(f,[lbl.dir 'stats_' lbl.analysis '.fig']);
close;
%}

end