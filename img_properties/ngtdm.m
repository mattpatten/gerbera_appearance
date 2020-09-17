function [ngtdm_stats] = ngtdm(imgs, param, lbls)

% The neighbourhood grey-tone difference matrix is a vector containing, for each gray level, a sum of the differences in 
% gray tone with all the surrounding pixels, for each pixel with that gray tone.
% 
% Algorithm parameters:
%     param.d       - The size of the neighbourhood to be considered. 
%                     1 indicates one pixel in each direction of the pixel in question (i.e., d=1 is 8 pixels in total, d=2 is 24 pixels).
%                     It is not recommended to use values greater than 2 as effects begin to wash out during averaging. 
%
% Inputs:
%     imgs.analysis - A series of masked, greyscale images with dimensions (x, y, 1(=colour), image_index)
%     imgs.mask     - For display purposes, lets us show the 'peripheral' pixels in one of the figures
%     param.d       - As above.
%     param.dispFig - Logical value for whether the function should display figure outputs.
%     lbls          - What the directory and filenames should be included to describe the specific analysis (.images, .colour, .analysis, .dir)
%
% Outputs:
%     ngtdm_stats   - The statistics generated from the NGTDM, as a struct. The following statistics are described (and quantified) in Amadasun & King (1989):
%          .coarseness - In a coarse texture, the primitives or basic patterns making up the texture are large, and the texture tends to possess a high degree 
%                        of local uniformity in intensity.
%          .contrast   - Perceptually, an image is said to have a high level of contrast if areas of different intensity levels are clearly visible. Thus high
%                        contrast means that the intensity difference between neighbouring regions is large.
%          .busyness   - A busy texture is one in which there are rapid changes of intensity from one pixel to its neighbour; that is the spatial frequency of 
%                        intensity changes is very high.
%          .complexity - This refers to the visual information content of a texture. A complex (or high information content) texture is one where there are many
%                        patches or primitives present in the texture, and more so when the primitives have different average intensities.
%          .strength   - A strong texture is one when the primitives that comprise it are easily definable and clearly visible. They tend to look attractive and 
%                        have a high degree of visual feel. 
% 
%
% Created by Matt Patten
% Created in September, 2018


%% Initialize

%constants
boxWidth = 2*param.d+1; %how many pixels wide we will be checking for the neighbourhood
greylevels = 0:255; %how many grey level intensities are we measuring.
[xSize, ySize, ~] = size(imgs.analysis);

ngtdm = NaN(size(imgs.analysis)); %pre-allocate


%% Compute NGTDM
ngtd_cols = im2col(imgs.analysis,[boxWidth boxWidth]); %load masked greyscale image
ctr = ceil(size(ngtd_cols,1)/2); %finds middle row of neighbourhood to be excluded (as this is the pixel itself)
ngtd_cols = ngtd_cols([1:(ctr-1) (ctr+1):end],:); %remove middle row (the pixel itself which should not be included in the neighbourhood),
ngtd_entries = round(1/(boxWidth^2-1) * sum(ngtd_cols,1,'includenan')); %sum each column and divide by the size of the neighbourhood
%put back in the shape of the original matrix, leaving space for the border (that can't be measured as it's not the full neighbourhood size)
ngtdm((1+param.d):(end-param.d),(1+param.d):(end-param.d)) = reshape(ngtd_entries, xSize-2*param.d,ySize-2*param.d);
this_ngtdm = ngtdm; %for indexing purposes later on
borders = this_ngtdm - imgs.analysis;
borders = borders / max(max(borders));

%generate imgs.mask from ngtdm (thus, excluding any edge regions - as they have been removed by including NaNs in the summation)
ngtdm_mask = double(~isnan(this_ngtdm)); %must convert to double to allow NaN's in place of 0's.
ngtdm_mask(ngtdm_mask==0)=NaN; %convert anything outside of the imgs.mask to zeros
ngtdm_pixels = round(imgs.analysis) .* ngtdm_mask; %apply imgs.mask to reduce image to exclude peripheral pixels (that are in neighbourhoods but whose values themselves are not counted)

%find indexes of each grey level
for i = 1:length(greylevels)
    this_grey = find(ngtdm_pixels==greylevels(i)); %find all the pixels of this specific intensity
    ngtdm_greycount(i) = length(this_grey); %how many of this pixel intensity is in the image
    
    if ngtdm_greycount(i) ~= 0
        ngtdm_tally(i) = sum(abs(i - this_ngtdm(this_grey))); %sum of the differences between grey value and for all pixels in the NGTDM (for all pixels of this grey)
    else %if there are no entries, set as zero
        ngtdm_tally(i) = 0;
    end
end

%% Compute statistics
%n2 = (xSize - 2*param.d) * (ySize-2*param.d); %number of pixels in matrix (non-periphery)
n2 = nansum(nansum(ngtdm_mask)); %number of (non periphery) pixels in matrix AFTER MASKING
Ng = sum(ngtdm_greycount>0); %how many different shades of grey are in the image (NOT the range)

%restrict greyscale range to only what exists in the current image
lims = [find(ngtdm_greycount,1,'first') find(ngtdm_greycount,1,'last')]; %compute greyscale range of image
greyRange = lims(1):lims(2); %scaling out extraneous grey values outside the range of those present in the picture
N = ngtdm_greycount(greyRange); %the range of grey values in the image
p = (N / n2); %proportion of pixels with this grey value over all pixels (excluding the periphery)
s = ngtdm_tally(greyRange); %the actual NGTDM: sum of the abs difference between pixel intensity and its neighbourhood for all pixels of a specific greyvalue

%coarseness
ngtdm_stats.coarseness = 1 / (sum(p .* s) + eps);

%contrast
sum_con = 0;
sum_con2 = 0;
for i=1:length(greyRange)
    for j=1:length(greyRange)
        sum_con = sum_con + p(i)*p(j)*(greyRange(i)-greyRange(j))^2;
    end
    %sum_con2 = sum_con2 + p .* circshift(p,i) .* (greyRange-circshift(greyRange,i)).^2; %same thing, different implementation. Just debugging to check.
end

ngtdm_stats.contrast = ( 1 / (Ng * (Ng - 1)) * sum_con ) * (1/n2 * sum(s));

%busyness
sum_busy = 0;
for i=1:length(greyRange)
    for j=1:length(greyRange)
        if p(i)~=0 && p(j)~=0 %neither of the grey tones can be empty
            %sum_busy = sum_busy + i*p(i) - j*p(j);
            sum_busy = sum_busy + abs(greyRange(i)*p(i) - greyRange(j)*p(j)); %absolute value isn't directly specified in eqn, but in text: 'magnitude of differences between...' and is used in example code
        end
    end
end

ngtdm_stats.busyness = sum(p .* s) / sum_busy;

%complexity
sum_complx = 0;
for i=1:length(greyRange)
    for j=1:length(greyRange)
        if p(i)~=0 && p(j)~=0 %neither of the grey tones can be empty
            sum_complx = sum_complx + (abs(greyRange(i) - greyRange(j)) / (n2 * (p(i) + p(j)))) * (p(i)*s(i) + p(j)*s(j));
        end
    end
end

ngtdm_stats.complexity = sum_complx;

%strength
sum_strength = 0;
for i=1:length(greyRange)
    for j=1:length(greyRange)
        if p(i)~=0 && p(j)~=0 %neither of the grey tones can be empty
            sum_strength = sum_strength + (p(i) + p(j)) * (greyRange(i) - greyRange(j))^2;
        end
    end
end

ngtdm_stats.strength = sum_strength / (eps + sum(s));

if param.dispFig
    %% Plotting for debugging purposes
    figure;
    set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
    %original image
    f = subplot(2,3,1);
    imshow(imgs.analysis/255);
    title('Original');
    axis square; axis off;
    
    %ngtdm (barely noticeable)
    f = subplot(2,3,2);
    imagesc(this_ngtdm);
    title('NGTDM (Scaled)');
    colormap('gray'); %colorbar;
    axis square; axis off;
    
    %Difference between ngtdm and the original image (scaled)
    f = subplot(2,3,3);
    imagesc(this_ngtdm - imgs.analysis);
    title('Diff (scaled): NGTDM - Orig');
    colormap('gray'); %colorbar;
    axis square; axis off;
    
    %'Peripheral' pixels excluded from primary use in NGTDM
    f = subplot(2,3,4);
    imagesc(imgs.mask - ~isnan(ngtdm_mask)); %original imgs.mask minus imgs.mask made from the indexes used in our analysis
    title('Peripheral pixels');
    colormap('gray');
    axis square; axis off;
    
    %ngtdm
    f = subplot(2,3,5);
    plot(0:255,ngtdm_greycount,'b','LineWidth',1);
    hold all;
    plot(0:255,ngtdm_tally,'r','LineWidth',1);
    title('Distribution of Greys');
    set(gca,'XLim',[0 255],'XTick',0:50:255);
    legend({'Num Greys','Tally (Grey - NGTD)'},'Location','NorthOutside');
    box off;
    
    %some tweaking of overall image properties
    set(gcf, 'Position', [1, 1, 900, 600]); %set size on screen
    saveas(f,[lbls.dir 'ngtdm_' lbls.image '_' lbls.analysis '.png']);
    close;
end
%{
    %% ATTEMPT AT USING FOR SEGMENTATION PURPOSES
    
    %binarize edges for edge detection/segmentation
    bw_borders = imbinarize(borders,'global'); %definitely global - adaptive is terrible
    red = cat(3,repmat(0.75,xSize,ySize),zeros(xSize,ySize),zeros(xSize,ySize));
    
    f = subplot(1,3,1);
    imshow(bw_borders); %picture of border alone
    title('NGTDM - Orig');
    
    f = subplot(1,3,2); hold all;
    imshow(imgs.analysis/255); %display in subplot
    imagesc(red,'AlphaData',bw_borders);
    title('Overlay');
    axis image; axis off;
    
    f = subplot(1,3,3);
    imshow(imgs.analysis/255);
    title('Original');
    
    %some tweaking of overall image properties
    set(gcf, 'Position', [1, 1, 1200, 500]); %set size on screen
    saveas(f,[lbls.dir 'segment_' lbls.image '_' lbls.analysis '.png']);
    close;
%}

%{
%% Display results
if param.dispFig
    
    %Display results - statistics
    figure;
    f = subplot(3,2,1);
    plot(ngtdm_stats.coarseness,'LineStyle','-','LineWidth',1,'Color','k');
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages);
    title('Coarseness');
    
    f = subplot(3,2,2);
    plot(ngtdm_stats.contrast,'LineStyle','-','LineWidth',1,'Color','k');
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages);
    title('Contrast');
    
    f = subplot(3,2,3);
    plot(ngtdm_stats.busyness,'LineStyle','-','LineWidth',1,'Color','k');
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages);
    title('Busyness');
    
    f = subplot(3,2,4);
    plot(ngtdm_stats.complexity,'LineStyle','-','LineWidth',1,'Color','k');
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages);
    title('Complexity');
    xlabel('Image #');
    
    f = subplot(3,2,5);
    plot(ngtdm_stats.strength,'LineStyle','-','LineWidth',1,'Color','k');
    set(gca,'XLim',[1 nImages],'XTick',1:10:nImages);
    title('Strength');
    xlabel('Image #');
    
    %some tweaking of overall image properties
    set(gcf, 'Position', [1, 1, 1000, 800]); %set size on screen
    saveas(f,[lbls.dir 'ngtdm_statistics_' lbls.analysis '.png']);
    saveas(f,[lbls.dir 'ngtdm_statistics_' lbls.analysis '.fig']);
    close;
end
%}

end