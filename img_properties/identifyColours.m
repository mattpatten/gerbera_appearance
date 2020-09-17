function [primaryLAB, secondaryLAB] = identifyColours(imgAbbr, gerb, segment_name, label, within_bounds, threshold)

% This generates colour swatch figures for the primary (and if applicable) secondary colours of the flower.
% It displays the RHS and UPOV codes and English descriptor of the colour(s). In addition, it generates histograms
% for the RHS and UPOV indexes.
%
% Inputs:
%    imgAbbr        - An abbrevation linking to the dataset of flowers to process.
%    gerb           - Image matrices for an individual flower, requiring specifically RGB, LAB, and Mask properties.
%    segment_name   - String to identify the segment of the flower we are processing here
%    label          - String providing a name/label for the flower.
%    within_bounds  - Flag to identify if flower likely contains two colours.
%    threshold      - the colour difference threshold that has to be exceeded for it to be considered to have a gradient
%
% Outputs:
%    primaryLAB   - The identified main L*a*b* colour of the ray florets for the flower
%    secondaryLAB - The identified secondary L*a*b* colour (if applicable) of the ray florets of the flower, NaN otherwise
%
%    In addition, figures containing a summary of the colours present in the petals of each flower.
%
% Created by Matt Patten in Jan, 2019


%get directories
[outputDir, imgPropDir] = get_dir(imgAbbr,'output','imgProp');

load([imgPropDir 'RHS_UPOV_Lab_colours.mat']); %load official RHS (royal horticultural society) colour labels in L*a*b* colour spectrum

%get properties
[xSize, ySize, ~] = size(gerb.RGB);

%generalize across segments into a single value
gerb.mask.seg1 = eval(['gerb.mask.' segment_name 'Seg1']);
gerb.mask.seg2 = eval(['gerb.mask.' segment_name 'Seg2']);
gerb.mask.seg  = eval(['gerb.mask.' segment_name]);
gerb.LAB_masked.seg  = eval(['gerb.LAB_masked.' segment_name]);

%replace anything outside the mask with NaNs
gerb.LAB_masked.seg1 = applyMask(gerb.LAB, gerb.mask.seg1);
gerb.LAB_masked.seg2 = applyMask(gerb.LAB, gerb.mask.seg2);

%if this flower has a gradient (i.e., consists of two separate colours)
%if within_bounds
    
    %draw main colour for each segment of petals
    [RHSvals1, RHSidx1, UPOVvals1, UPOVidx1] = computeMainColour(gerb.LAB_masked.seg1, lab_values, UPOV_labels);
    [RHSvals2, RHSidx2, UPOVvals2, UPOVidx2] = computeMainColour(gerb.LAB_masked.seg2, lab_values, UPOV_labels);
    
    %swap values of variables so whichever segment has the mode with a higher frequency is set as the first one
    if RHSvals2(1) >= RHSvals1(1)
        [RHSvals2,RHSvals1]   = deal(RHSvals1,RHSvals2); %swaps values of variables
        [RHSidx2,RHSidx1]     = deal(RHSidx1,RHSidx2);
        [UPOVvals2,UPOVvals1] = deal(UPOVvals1,UPOVvals2);
        [UPOVidx2,UPOVidx1]   = deal(UPOVidx1,UPOVidx2);
        [gerb.mask.seg2,gerb.mask.seg1] = deal(gerb.mask.seg1,gerb.mask.seg2);
    end
    
    %defined for readability
    primaryRHS = RHSidx1(1);
    primaryUPOV = UPOVidx1(1);
    
    %find first secondary colour that is a different UPOV label than from the primary
    if UPOVidx2(1)~=primaryUPOV
        secondaryUPOV = UPOVidx2(1);
    elseif UPOVidx2(1)==primaryUPOV
        secondaryUPOV = UPOVidx2(2);
    end
    
    %find all RHS values that match this UPOV colour, and find the match with the smallest index / closest to the top / most frequent
    [~,match_idx] = ismember(find(UPOV_labels==secondaryUPOV),RHSidx2);
    secondaryRHS = RHSidx2(min(match_idx));
    
    %primary and secondary colours - create colour swatches (from mode)
    colourpatch1 = repmat(reshape(lab2rgb(lab_values(primaryRHS,:)),1,1,3),xSize,ySize);
    colourpatch2 = repmat(reshape(lab2rgb(lab_values(secondaryRHS,:)),1,1,3),xSize,ySize);
    
    %combine segments for histogram
    UPOVvals1(UPOVidx1,1) = UPOVvals1;
    UPOVvals2(UPOVidx2,1) = UPOVvals2;
    UPOVvals = UPOVvals1 + UPOVvals2;
    RHSvals1(RHSidx1,1) = RHSvals1;
    RHSvals2(RHSidx2,1) = RHSvals2;
    RHSvals = RHSvals1 + RHSvals2;
    
    %save for output
    primaryLAB   = lab_values(primaryRHS,:);
    secondaryLAB = lab_values(secondaryRHS,:);
    
%else %treat petals as a single segment
%    [RHSvals, RHSidx, UPOVvals, UPOVidx] = computeMainColour(gerb.LAB_masked.seg, lab_values, UPOV_labels);
%    primaryRHS = RHSidx(1);
%    colourpatch1 = repmat(reshape(lab2rgb(lab_values(primaryRHS,:)),1,1,3),xSize,ySize);
%    
%    %re-order for histogram
%    UPOVvals(UPOVidx,1) = UPOVvals;
%    RHSvals(RHSidx,1) = RHSvals;
%    
%    %save for output
%    primaryLAB   = lab_values(primaryRHS,:);
%    secondaryLAB = lab_values(primaryRHS,:);
%end

%sort where pictures are saved depending on threshold
colourDiff = pdist([primaryLAB; secondaryLAB]);
if and(colourDiff > threshold, within_bounds)
    colourPatchDir = [outputDir filesep 'colourPatch_' segment_name filesep]; 
    if ~exist(colourPatchDir,'dir'), mkdir(colourPatchDir); end %output directory for this analysis
else
    colourPatchDir = [outputDir filesep 'colourPatch_' segment_name '_excluded' filesep]; 
    if ~exist(colourPatchDir,'dir'), mkdir(colourPatchDir); end %output directory for this analysis
end

f = figure; nRows = 3; nCols = 3;
set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things

%plot original image
subplot(nRows,nCols,1);
imshow(gerb.RGB);
title(['Diff: ' num2str(colourDiff)]);
%title('Original');

%if within_bounds
    %draw figure for segment of petals
    subplot(nRows,nCols,2); imshow(applyMask(gerb.RGB, gerb.mask.seg1));
    subplot(nRows,nCols,3); imshow(applyMask(gerb.RGB, gerb.mask.seg2)); 
    
    %draw colour patches of the 2 most frequent pixel colours
    subplot(nRows,nCols,5); imshow(colourpatch1); title(['RHS: ' RHS_labels{primaryRHS}   ',  UPOV: ' num2str(UPOV_labels(primaryRHS))]);   xlabel(colour_labels{primaryRHS});
    subplot(nRows,nCols,6); imshow(colourpatch2); title(['RHS: ' RHS_labels{secondaryRHS} ',  UPOV: ' num2str(UPOV_labels(secondaryRHS))]); xlabel(colour_labels{secondaryRHS});
    
%else
%    %draw petals and calculated colour
%    subplot(nRows,nCols,2); imshow(applyMask(gerb.RGB, gerb.mask.seg));
%    subplot(nRows,nCols,5); imshow(colourpatch1); title(['RHS: ' RHS_labels{primaryRHS} ',  UPOV: ' num2str(UPOV_labels(primaryRHS))]); xlabel(colour_labels{primaryRHS});
%end

%UPOV and RHS histograms, respectively
subplot(nRows,nCols,4); bar(1:length(UPOVvals),UPOVvals); title('Histogram'); xlim([0 size(UPOVvals,1)]); xlabel('UPOV Labels Idx'); xticks(0:10:length(UPOVvals)); ylabel('Freq'); box off;
subplot(nRows,nCols,7); bar(1:length(RHSvals), RHSvals);  title('Histogram'); xlim([0 size(RHSvals,1)]);  xlabel('RHS Labels Idx');  ylabel('Freq'); box off;

set(gcf, 'Position', [1, 1, 1200, 1000]); %set size on screen
saveas(f,[colourPatchDir 'colourPatch_' label '.png']);
close;

end


function [RHSfreq, RHSidx, UPOVfreq, UPOVidx] = computeMainColour(data, lab_values, UPOV_values)

%flatten each channel of L*a*b* into a column vector (petals only)
flatlab = reshape(data,numel(data)/size(data,3),size(data,3)); %each channel becomes a column
flatlab = flatlab(any(~isnan(flatlab),2),:); %removes nans to reduce array size

%get Euclidean distance of each pixel to all RHS-defined colours and choose whichever 'colour' is closest
closestRHS = zeros(size(flatlab,1),1);
for i=1:size(flatlab,1)
    [~, closestRHS(i,:)] = min(vecnorm(flatlab(i,:) - lab_values,2,2));
end

%find RHS-defined colour that is closest to the most amount of pixels
[RHSfreq, RHSidx] = maxk(histcounts(closestRHS,size(lab_values,1),'BinLimits',[1 size(lab_values,1)]),size(lab_values,1));
RHSfreq = RHSfreq'; RHSidx = RHSidx'; %transpose to column vector

%combine different RHS values that link to the same UPOV label
UPOVfreq = arrayfun(@(x) sum(RHSfreq(UPOV_values(RHSidx)==x)),unique(UPOV_values)); %for each UPOV value (1:73), sum the freq of all RHS matches for this UPOV
[UPOVfreq, UPOVidx] = sort(UPOVfreq,'descend');

end
