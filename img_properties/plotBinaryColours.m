function [params] = plotBinaryColours(imgAbbr, gerb, segment_name, has_gradient, primaryColour, secondaryColour, seg_start, seg_end, label, limits)

% This takes the primary and secondary colours of any flowers identified to have a gradient and re-plots the ray florets using these two colours only.
% It then fits a psychometric function to the rate of change between the two colours (as a function of eccentricity) to identify the point where the
% gradient changes (threshold), and whether it changes rapidly or gradually (slope).
%
% Inputs:
%    imgAbbr          - Dataset abbreviation title.
%    gerb             - A struct of all the different formats for this particular flower.
%    segment_name     - String to identify the segment of the flower we are processing here
%    has_gradient     - Boolean identifying whether current flower has a gradient.
%    primaryColour    - The L*a*b* values of the primary colour of the ray florets for this flower.
%    secondaryColour  - The L*a*b* values of the secondary colour of the ray florets for this flower.
%    seg_start        - The eccentricity where this segment of the flower starts (closest to the centre).
%    seg_end          - The eccentricity where this segment of the flower ends (closest to the outside).
%    segmentation     - The eccentricities separating different flower sections (disk to trans, trans to ray, ray to the outside edge).
%    label            - String providing a name/label for this flower.
%    limits           - Bounds that the threshold must be within for it to be considered an actual change rather than the edge of the segment and the start of the next segment, e.g., [0.05 0.95] 
%
% Outputs:
%    params           - A cell array containing 4 values: the threshold, the slope and the min/max of the above described psychometric function.
%
% Created by Matt Patten in Feb, 2019


%File I/O
outputDir = get_dir(imgAbbr,'output');



%get properties
[xSize, ySize, nColours] = size(gerb.RGB);
%params = cell(nImages,1);
%params = arrayfun(@(~) zeros(1,4), (1:nImages)', 'UniformOutput', false); %create [0,0,0,0] for each cell, one for each flower
params = zeros(1,4); %create [0,0,0,0]
options = optimset('Display','none');

mask_segment = eval(['gerb.mask.' segment_name]);
    
if has_gradient
    %initialize
    maskVals = [];
    primImg = zeros(xSize,ySize);
    secImg  = zeros(xSize,ySize);
    
    %get indexes of values inside the petals
    maskIdxs = find(mask_segment==1);
    
    %add each colour channel as a column vector (for indexes inside the mask only)
    for cc=1:nColours
        tmp = gerb.LAB(:,:,cc);
        maskVals = [maskVals tmp(maskIdxs)];
    end
    
    %identify which of the two colours is closer to each pixel
    [~, minClass] = min([vecnorm(maskVals - primaryColour,2,2) vecnorm(maskVals - secondaryColour,2,2)],[],2);
    minClass = 2 - minClass; %sets primary colour to one, secondary colour to zero
    
    %Split the indexes of the mask into which pixels preferred which colour
    primIdxs = maskIdxs(minClass==1);
    secIdxs  = maskIdxs(minClass==0);
    
    %get the indexes of the original mask (so we can feed it back into the original picture at the right location)
    primImg(primIdxs)=1;
    secImg(secIdxs)=1;
    
    %change from binary (0,0,0) and (1,1,1) to the primary and secondary colour values
    primColImg = cat(3, primImg * primaryColour(1), primImg * primaryColour(2), primImg * primaryColour(3)); %just the primary colour
    secColImg = cat(3, secImg * secondaryColour(1), secImg * secondaryColour(2), secImg * secondaryColour(3)); %just the secondary colour
    bothColImg = primColImg + secColImg; %overlay both colours onto the one image
    
    %get the (stable) fminsearchbnd parameters
    xvals = (seg_start+1):seg_end;
    fit_steps = (seg_start+1):0.1:seg_end; %smaller steps for smoother line
    init_params = [(seg_end + seg_start) / 2, 1, 0, 1]; %threshold, slope, min, max
    lower_bound = [(seg_start+1) 0 0 0]; %threshold, slope, min, max
    upper_bound = [seg_end Inf 1 1];   %threshold, slope, min, max
    
    %get proportion of pixels per eccentricity (and inside the mask) belonging to the primary colour
    masked_flower = applyMask(primImg, mask_segment); %NaNs outside of mask, ones only for primary colour
    propPrim = findRadialAvgs(masked_flower,1,round(xSize/2),round(ySize/2),'sum',1); %get proportion of ones in each eccentricity
    propPrim = propPrim((seg_start+1):seg_end); %reduce array size to just ray florets
    paramsPrim = fminsearchbnd(@(paramsPrim) compute_fit_residuals(paramsPrim,xvals,propPrim),init_params,lower_bound,upper_bound,options); %find psychometric/logistic parameters that minimize the residuals
    
    %get proportion of pixels per eccentricity (and inside the mask) belonging to the secondary colour
    masked_flower = applyMask(secImg, mask_segment); %NaNs outside of mask, ones only for secondary colour
    propSec = findRadialAvgs(masked_flower,1,round(xSize/2),round(ySize/2),'sum',1); %get proportion of ones in each eccentricity
    propSec = propSec((seg_start+1):seg_end); %reduce array size to just ray florets
    paramsSec = fminsearchbnd(@(paramsSec) compute_fit_residuals(paramsSec,xvals,propSec),init_params,lower_bound,upper_bound,options); %find psychometric/logistic parameters that minimize the residuals
    
    if compute_fit_residuals(paramsPrim,xvals,propPrim) < compute_fit_residuals(paramsSec,xvals,propSec) %whichever is a better fit to the data
        fit_prop_prim = logistic4(paramsPrim,fit_steps); %fit optimal parameters to function
        fit_prop_sec  = 1 - fit_prop_prim; %other function is inverse of this
        
        %regardless of which direction the gradient change is (captured by the primary and secondary colours), positive values indicate the structure of the gradient change
        params = [paramsPrim(1) abs(paramsPrim(2)) paramsPrim(3) paramsPrim(4)];
        
    else
        fit_prop_sec  = logistic4(paramsSec,fit_steps); %fit optimal parameters to function
        fit_prop_prim = 1 - fit_prop_sec; %other function is inverse of this
        
        %regardless of which direction the gradient change is (captured by the primary and secondary colours), positive values indicate the structure of the gradient change
        params = [paramsSec(1) abs(paramsSec(2)) paramsSec(3) paramsSec(4)];
    end
    
    [~, xThrIdx] = min(abs(params(1) - fit_steps)); %find index of line where threshold is (for plotting vertical line of threshold)
    
    %normalize threshold value (if nonzero)
    unnormed_threshold = params(1);
    if params(1)~=0 
        params(1) = (params(1) - seg_start)/(seg_end - seg_start);
    end
    
    %draw figures
    f = figure; rows = 2; cols = 5;
    set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
    
    subplot(rows,cols,1);
    imshow(gerb.RGB);
    title('Original');
    
    subplot(rows,cols,2);
    imshow(lab2rgb(primColImg));
    title('Primary Colour'); 
    xlabel(['Threshold: ' num2str(params(1))]);
    
    subplot(rows,cols,3);
    imshow(lab2rgb(secColImg));
    title('Secondary Colour');
    
    subplot(rows,cols,4);
    imshow(lab2rgb(bothColImg));
    title('Both');
    xlabel(['Diff: ' num2str(pdist([primaryColour; secondaryColour]))]);
    
    subplot(rows,cols,5);
    imshow(applyMask(gerb.RGB, mask_segment));
    title([segment_name ' florets']);
    
    %psychometric function - primary colour
    subplot(rows,cols,7); hold all;
    plot([unnormed_threshold unnormed_threshold],[0 fit_prop_prim(xThrIdx)],'Marker','none','LineStyle','--','LineWidth',1,'Color',rgb('Grey')); %threshold line
    plot(fit_steps,fit_prop_prim,'Marker','none','LineStyle','-','LineWidth',1,'Color',rgb('Red')); %psychometric fit
    plot(xvals,propPrim,'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor',rgb('Blue'),'MarkerSize',3); %data values
    axis square; box off; xlim([seg_start seg_end]); ylim([0 1]);
    xlabel('Eccentricity'); ylabel('Number of pixels (prop)');
    title('Primary Colour Pixels');
    
    %psychometric function - secondary colour
    subplot(rows,cols,8); hold all;
    plot([unnormed_threshold unnormed_threshold],[0 fit_prop_sec(xThrIdx)],'Marker','none','LineStyle','--','LineWidth',1,'Color',rgb('Grey')); %threshold line
    plot(fit_steps,fit_prop_sec,'Marker','none','LineStyle','-','LineWidth',1,'Color',rgb('Red')); %psychometric fit
    plot(xvals,propSec,'LineStyle','none','Marker','o','MarkerEdgeColor','None','MarkerFaceColor',rgb('Blue'),'MarkerSize',3); %data values
    axis square; box off; xlim([seg_start seg_end]); ylim([0 1]);
    xlabel('Eccentricity'); ylabel('Number of pixels (prop)');
    title('Secondary Colour Pixels');
    
    %only include if not too close to the edge (probably is edge of next segment) - otherwise, reset
    if and(params(1)>limits(1),params(1)<limits(2))
        binaryColsDir = [outputDir 'binarized_colours_' segment_name filesep];
        if ~exist(binaryColsDir,'dir'), mkdir(binaryColsDir); end

    elseif or(params(1)<limits(1),params(1)>limits(2))
        binaryColsDir = [outputDir 'binarized_colours_' segment_name '_excluded' filesep];
        if ~exist(binaryColsDir,'dir'), mkdir(binaryColsDir); end
        params = zeros(1,4); %reset to [0,0,0,0]
    end
        
    set(gcf, 'Position', [1, 1, 1600, 600]); %set size on screen
    saveas(f,[binaryColsDir 'binaryCols_' label '.png']);
    close;

end

end