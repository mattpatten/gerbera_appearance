function get_flower_properties(imgAbbr)

% Extracts numerical details of the flower for future analysis
%
% Created by Matt Patten
% Created in January, 2019

clear global;
global feature_count;
global imgProperties;
global header;

%% Set up file I/O
disp('Loading images....');

%get/create directories
[imgDir, maskDir, outputDir, dataDir] = get_dir(imgAbbr,'img','mask','output','data');

%load image details from labels file
load([imgDir 'labels_' imgAbbr '.mat'],'labels');
try load([imgDir 'segmentation_manual.mat'],'segmentation'); %load manual segmentations (if file is there)
catch, load([imgDir 'segmentation.mat'],'segmentation');     %otherwise load automatically-generated version
end

% Get flwrpoll survey responses
%read_in_flowerpoll_attributes(imgAbbr); %generates mat files in data directory with relevant properties
load([dataDir 'flwrpoll_ratings_' imgAbbr '_physical.mat'],'poll');


fprintf('\n'); %new line for user updates
for flowerIdx = 1:length(labels)

    fprintf('Getting flower properties for flower #%d: ', flowerIdx); %update user
    
    %load image and mask for this flower
    gerb.RGB        = imread([imgDir  labels{flowerIdx}  '.png']); %load image
    gerb.mask.whole = imread([maskDir labels{flowerIdx} 'm.png']); %load mask
    seg = segmentation(flowerIdx,:); %load segmentation values

    imageID = str2num(regexp(labels{flowerIdx},'\d.*','match','once')); %get integer related to image filename label
    
    %% Get some basic image properties
    [xSize, ySize, nColours] = size(gerb.RGB);
    nSegments = 3;
    r = convert2polar([xSize ySize]); 
    imRad = floor(xSize/2);
    xCtr  = round(xSize/2); %col
    yCtr  = round(ySize/2); %row
    
    
    %% Divide flower into segments
    gerb.mask.disk  = r<=seg(1); %centre / disk florets
    gerb.mask.trans = and(r>seg(1),r<=seg(2)); %middle / trans florets
    gerb.mask.ray   = and(r>seg(2),gerb.mask.whole); %petals / ray florets
    
    %compute ray widths / petal lengths
    outsideMask = findRadialAvgs(~gerb.mask.whole,1,round(xSize/2),round(ySize/2),'sum',0); %sum anything that's outside of mask (per eccentricity)
    seg(3) = find(outsideMask>0.5,1); %add outside edge eccentricity as third "segmentation"
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Performing transforms      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Transforming / ');
    
    %% Colour transforms
    %gerb.HSV = rgb2hsv(gerb.RGB);
    gerb.LAB  = rgb2lab(gerb.RGB);
    gerb.grey = rgb2gray(gerb.RGB);
    
    %{
    %% Polar transforms (and polar mask segmentation)
    gerb.polGrey.whole = polartrans(gerb.grey,imRad,360,xCtr,yCtr,'linear','valid'); %greyscale
    gerb.polGrey.ray = gerb.polGrey.whole;
    gerb.polGrey.ray(1:seg(2),:) = 0; %cut mask down to just ray florets
    gerb.polMask.whole = imbinarize(polartrans(gerb.mask.whole,imRad,360,xCtr,yCtr,'linear','valid')); %mask
    
    for cc=1:3 %colour transforms
        gerb.polRGB(:,:,cc) = polartrans(gerb.RGB(:,:,cc),imRad,360,xCtr,yCtr,'linear','valid'); %RGB file
        gerb.polLAB(:,:,cc) = polartrans(gerb.LAB(:,:,cc),imRad,360,xCtr,yCtr,'linear','valid'); %LAB file
    end
    
    %polar mask segmentation for ray florets
    gerb.polMask.ray = gerb.polMask.whole;
    gerb.polMask.ray(1:(seg(2)-1),:) = 0;
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Applying masks to transforms   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Applying mask / ');

    %apply masks to transformations
    gerb.LAB_masked.whole = applyMask(gerb.LAB, gerb.mask.whole); %LAB
    gerb.LAB_masked.disk  = applyMask(gerb.LAB, gerb.mask.disk);
    gerb.LAB_masked.trans = applyMask(gerb.LAB, gerb.mask.trans);
    gerb.LAB_masked.ray   = applyMask(gerb.LAB, gerb.mask.ray);
    
    gerb.grey_masked.whole = applyMask(double(gerb.grey), gerb.mask.whole); %Grey
    
    %{
    gerb.polLAB_masked.ray  = applyMask(gerb.polLAB, gerb.polMask.ray); %Polar
    gerb.polRGB_masked.ray  = applyMask(gerb.polRGB, gerb.polMask.ray);
    gerb.polGrey_masked.ray = applyMask(gerb.polGrey.whole,gerb.polMask.ray);
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Descriptive statistics on polar transform   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    fprintf('Descriptives / ');
    
    img.polarDescr.analysis = gerb.polLAB_masked.ray; %compute polar descriptive statistics on ray florets only
    img.polarDescr.display  = gerb.RGB;
    plotPolarDecriptives(img.polarDescr,imgAbbr, labels{flowerIdx});
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Determining petal edge circularity   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    fprintf('Edge circularity / ');

    img.ang.analysis       = gerb.polGrey.ray;
    img.ang.mask           = gerb.polMask.whole;
    img.ang.regularDisplay = gerb.RGB;
    img.ang.polarDisplay   = gerb.polRGB;

    getCircularity(img.ang, imgAbbr, labels{flowerIdx});
    %}
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Radial Analysis         %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Radial analysis / ');
    
    %compute mean radial values
    radmean_LAB.whole(:,1) = findRadialAvgs(gerb.LAB_masked.whole(:,:,1),1,round(xSize/2),round(ySize/2),'mean',1); %get mean value for each colour channel at each eccentricity
    radmean_LAB.whole(:,2) = findRadialAvgs(gerb.LAB_masked.whole(:,:,2),1,round(xSize/2),round(ySize/2),'mean',1);
    radmean_LAB.whole(:,3) = findRadialAvgs(gerb.LAB_masked.whole(:,:,3),1,round(xSize/2),round(ySize/2),'mean',1);
    
    %replicate all eccentricities for each of the segments
    radmean_LAB.disk  = radmean_LAB.whole;
    radmean_LAB.trans = radmean_LAB.whole;
    radmean_LAB.ray   = radmean_LAB.whole;
    
    %limit analysis to specific segments (i.e., remove values for eccentricities outside of ray florets)
    radmean_LAB.disk (            (seg(1)+1):end ,:) = nan;
    radmean_LAB.trans([1:(seg(1)) (seg(2)+1):end],:) = nan;
    radmean_LAB.ray  ([1:(seg(2)) (seg(3)+1):end],:) = nan;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Segmentation of each segment  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Segmenting ray florets / ');

    %do not consider a change close to the edge of the segment as a proper gradient change (as it may be the start of the next area)
    edge_boundaries_ray   = [0.1 0.95];
    edge_boundaries_trans = [0.3 0.7];
    %edge_boundaries_disk  = [0   0.8];
    
   
    %separate ray florets into two parts
    [gerb.mask.raySeg1,   gerb.mask.raySeg2,   kmeans_idxs_ray]   = segmentFloretsRadially(radmean_LAB.ray, gerb.mask.ray); %adds gerb.mask.raySeg1 and gerb.mask.raySeg2 to gerb struct
    [gerb.mask.transSeg1, gerb.mask.transSeg2, kmeans_idxs_trans] = segmentFloretsRadially(radmean_LAB.trans, gerb.mask.trans); 
    %[gerb.mask.diskSeg1,  gerb.mask.diskSeg2,  kmeans_idxs_disk]  = segmentFloretsRadially(radmean_LAB.disk, gerb.mask.disk); 
        
    %if clustering cames back with more than one boundary, find segment with higher number of included eccentricities
    if length(kmeans_idxs_ray) > 1
        [~,idx] = max([kmeans_idxs_ray(1) - seg(2), seg(3) - kmeans_idxs_ray(2)]);
        ray_segment = kmeans_idxs_ray(idx);
    else
        ray_segment = kmeans_idxs_ray;
    end
    norm_ray_segment = (ray_segment - seg(2))/(seg(3) - seg(2)); %normalize threshold value (if nonzero)
    within_bounds_ray = and(norm_ray_segment > edge_boundaries_ray(1),norm_ray_segment < edge_boundaries_ray(2));
    
    if length(kmeans_idxs_trans) > 1
        [~,idx] = max([kmeans_idxs_trans(1) - seg(1), seg(2) - kmeans_idxs_trans(2)]);
        trans_segment = kmeans_idxs_trans(idx);
    else
        trans_segment = kmeans_idxs_trans;
    end
    norm_trans_segment = (trans_segment - seg(1))/(seg(2) - seg(1)); %normalize threshold value (if nonzero)
    within_bounds_trans = and(norm_trans_segment > edge_boundaries_trans(1), norm_trans_segment < edge_boundaries_trans(2));
    
    %{
    if length(kmeans_idxs_disk) > 1
        [~,idx] = max([kmeans_idxs_disk(1), seg(1) - kmeans_idxs_disk(2)]);
        disk_segment = kmeans_idxs_disk(idx);
    else
        disk_segment = kmeans_idxs_disk;
    end
    norm_disk_segment = disk_segment/seg(1); %normalize threshold value (if nonzero)
    within_bounds_disk = and(norm_disk_segment > edge_boundaries_disk(1), norm_disk_segment < edge_boundaries_disk(2));
    %}
    
    %limit analysis to each segment of the ray florets (and remove all other values)
    radmean_LAB.raySeg1   = radmean_LAB.whole;
    radmean_LAB.raySeg2   = radmean_LAB.whole;
    radmean_LAB.transSeg1 = radmean_LAB.whole;
    radmean_LAB.transSeg2 = radmean_LAB.whole;
    %radmean_LAB.diskSeg1  = radmean_LAB.whole;
    %radmean_LAB.diskSeg2  = radmean_LAB.whole;
    
    radmean_LAB.raySeg1([1:(seg(2)) (ray_segment+1):end],:)     = nan;
    radmean_LAB.raySeg2([1:ray_segment  (seg(3)+1):end],:)      = nan;

    radmean_LAB.transSeg1([1:(seg(1)) (trans_segment+1):end],:) = nan;
    radmean_LAB.transSeg2([1:trans_segment  (seg(2)+1):end],:)  = nan;
    
    %radmean_LAB.diskSeg1((disk_segment+1):end,:)                = nan;
    %radmean_LAB.diskSeg2([1:disk_segment  (seg(1)+1):end],:)    = nan;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Correlation analysis      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Correlation / ');
    
    %pre-allocate cell array
    xCorrVals = cell(5,1);
    yCorrFit  = cell(nColours,5); 
    coeffsLAB = cell(nColours,5); 
    corrR     = cell(nColours,5); 

    %correlation analysis (for slope of ray florets)
    for cc=1:nColours
        [xCorrVals{1}, yCorrFit{cc,1}, coeffsLAB{cc,1}, corrR{cc,1}] = correlationAnalysis(radmean_LAB.disk(:,cc));    %disk
        [xCorrVals{2}, yCorrFit{cc,2}, coeffsLAB{cc,2}, corrR{cc,2}] = correlationAnalysis(radmean_LAB.trans(:,cc));   %trans
        [xCorrVals{3}, yCorrFit{cc,3}, coeffsLAB{cc,3}, corrR{cc,3}] = correlationAnalysis(radmean_LAB.ray(:,cc));     %ray - all
    end
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Colour identification      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Gradients / ');
    
    colour_separation_threshold_ray   = 23;
    colour_separation_threshold_trans = 35;
    %colour_separation_threshold_disk  = 40;
    
    %compute the difference in colours between the two segments
    %ray_colourDist   = get_colour_distance(gerb.LAB, gerb.mask.raySeg1, gerb.mask.raySeg2);
    %trans_colourDist = get_colour_distance(gerb.LAB, gerb.mask.transSeg1, gerb.mask.transSeg2);
    %disk_colourDist  = get_colour_distance(gerb.LAB, gerb.mask.diskSeg1, gerb.mask.diskSeg2);
    
    %determining primary and (sometimes) secondary colour of each segment
    %[primaryColour, secondaryColour] = identifyColours(imgAbbr, gerb, labels{flowerIdx}, rays_grad_bool, 1, colour_separation_threshold_ray);
    [ray_primaryColour, ray_secondaryColour]     = identifyColours(imgAbbr, gerb, 'ray', labels{flowerIdx}, within_bounds_ray, colour_separation_threshold_ray);
    [trans_primaryColour, trans_secondaryColour] = identifyColours(imgAbbr, gerb, 'trans', labels{flowerIdx}, within_bounds_trans, colour_separation_threshold_trans);
    %[disk_primaryColour, disk_secondaryColour]   = identifyColours(imgAbbr, gerb, 'disk', labels{flowerIdx}, within_bounds_disk, colour_separation_threshold_disk);
    
    %find class of data with largest Euclidean distance between the primary and secondary colours
    ray_colourDiff = pdist([ray_primaryColour; ray_secondaryColour]);
    ray_grad_bool = and(ray_colourDiff > colour_separation_threshold_ray, within_bounds_ray);
    
    trans_colourDiff = pdist([trans_primaryColour; trans_secondaryColour]);
    trans_grad_bool = and(trans_colourDiff > colour_separation_threshold_trans, within_bounds_trans);
    
    %disk_colourDiff = pdist([disk_primaryColour; disk_secondaryColour]);
    %disk_grad_bool = and(disk_colourDiff > colour_separation_threshold_disk, within_bounds_disk);
    
    %convert flowers with gradients into two-colour binary image (otherwise return zeroes)
    ray_binPsychoParams   = plotBinaryColours(imgAbbr, gerb, 'ray',     ray_grad_bool,   ray_primaryColour,   ray_secondaryColour, seg(2), seg(3), labels{flowerIdx}, edge_boundaries_ray);
    trans_binPsychoParams = plotBinaryColours(imgAbbr, gerb, 'trans', trans_grad_bool, trans_primaryColour, trans_secondaryColour, seg(1), seg(2), labels{flowerIdx}, edge_boundaries_trans);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%       Generate L*a*b* histograms       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Histogram / ');
    
    plot_histograms_LAB(imgAbbr, gerb, seg, radmean_LAB.whole, [kmeans_idxs_ray kmeans_idxs_trans], xCorrVals, yCorrFit, coeffsLAB, corrR, labels{flowerIdx});

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     Grey-level Co-occurrence Matrix      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('GLCM / ');
    
    %define the set of images we are going to use for this analysis
    img.glcm.analysis = gerb.grey_masked.whole;
    img.glcm.display = gerb.RGB;
    
    %define parameters
    p.glcm.bins = 64; %the number of grey levels to consider (and consequently the size of the GLCM)
    p.glcm.symmetry = true; %whether we consider an inverse relationship. i.e., [0,1] as well as [1,0]
    p.glcm.grayLimits = []; %[0 255]; %image instensity values to scale between (i.e., min and max for greyscale values to consider, numbins will be divided within these values)
    p.glcm.offsetValues = 1:50; %distance to measure between pixels - can display results of multiple offsets
    
    %generate file structure
    names.glcm.analysis = ['GLCM_' num2str(p.glcm.bins) 'bins'];
    names.glcm.dir = [outputDir names.glcm.analysis filesep];
    names.glcm.image = labels{flowerIdx};
    if ~exist(names.glcm.dir,'dir'), mkdir(names.glcm.dir); end

    %Calculate Tamura features
    TamuraStats = Tamura3Sigs(img.glcm.analysis); %Coarseness / Coarseness_hist (1x3) / Directionality / Contrast
    %fprintf('Flower #%i (Coarseness, Directionality, Contrast): %.3f\t%.3f\t%.3f \n',flowerIdx,TamuraStats(flowerIdx,1),TamuraStats(flowerIdx,5),TamuraStats(flowerIdx,6));
    
    %Calculate Haralick features (functional - but doesn't have clear psychophysical comparison)
    %HaralickStats = compute_Haralick_statistics(img.glcm, p.glcm, names.glcm.analysis, 0);
    
    %An earlier version that displays GLCMs for various offsets - not used for calculations, but can be useful visually.
    glcm_stats = greyLevelCooccurrence(img.glcm, p.glcm, names.glcm); %display GLCMs
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                  NGTDM                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('NGTDM / ');
    
    
    %define the set of images we are going to use for this analysis
    img.ngtdm.analysis = double(gerb.grey_masked.whole); %gerbsLABMasked(:,:,3,:);
    img.ngtdm.mask = gerb.mask.whole;
    
    %define parameters
    p.ngtdm.d = 1; %size of neighbourhood (how many pixels wide will we check). Recommended 1 or 2 only.
    p.ngtdm.dispFig = 1; %whether to display figures with output
    
    %generate file structure
    names.ngtdm.image = labels{flowerIdx};
    names.ngtdm.analysis = ['v4_d' num2str(p.ngtdm.d)]; %label to use for save directory and in filename
    
    names.ngtdm.dir = [outputDir 'ngtdm_' names.ngtdm.analysis filesep];
    if ~exist(names.ngtdm.dir,'dir'), mkdir(names.ngtdm.dir); end %if directory doesn't exist, create it
    
    %perform analysis
    ngtdm_stats = ngtdm(img.ngtdm, p.ngtdm, names.ngtdm);
    
    %{
    %other version for comparison
    allgreys = unique(round(img.ngtdm.analysis));
    allgreys = allgreys(~isnan(allgreys));
    [otherNGTDM, validCount] = getNGTDM(img.ngtdm.analysis,allgreys);
    allstats = getNGTDMtextures(otherNGTDM,validCount);
    other_ngtdm_stats.coarseness = allstats.Coarseness;
    other_ngtdm_stats.contrast   = allstats.Contrast;
    other_ngtdm_stats.busyness   = allstats.Busyness;
    other_ngtdm_stats.complexity = allstats.Complexity;
    other_ngtdm_stats.strength   = allstats.Strength;
    %}
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%       Generate golden ratio info       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    gr = (1+sqrt(5))/2;
    gr_ss = (1 - seg(2)/(seg(1)*gr))^2 + (1 - seg(3)/(seg(2)*gr))^2;
    %gr_ss(gr_ss>10)=10; %implement maximum SS
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Extract survey responses         %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    poll_attributes = zeros(1,length(poll.headers)); %initialize
    
    [~, idx] = ismember(imageID,poll.flowerIDs); %find if, and where
    if idx>0 %if found
        poll_attributes = poll.means(:,idx);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Adding properties to table         %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    feature_count = 0;
    
    %% Colour of each segment
    add_feature(flowerIdx, 'colour_mean_disk_L',  nanmean(nanmean(gerb.LAB_masked.disk(:,:,1),1),2));
    add_feature(flowerIdx, 'colour_mean_disk_A',  nanmean(nanmean(gerb.LAB_masked.disk(:,:,2),1),2));
    add_feature(flowerIdx, 'colour_mean_disk_B',  nanmean(nanmean(gerb.LAB_masked.disk(:,:,3),1),2));
    
    add_feature(flowerIdx, 'colour_primary_trans_L', trans_primaryColour(1));
    add_feature(flowerIdx, 'colour_primary_trans_A', trans_primaryColour(2));
    add_feature(flowerIdx, 'colour_primary_trans_B', trans_primaryColour(3));

    add_feature(flowerIdx, 'colour_primary_ray_L',  ray_primaryColour(1));
    add_feature(flowerIdx, 'colour_primary_ray_A',  ray_primaryColour(2));
    add_feature(flowerIdx, 'colour_primary_ray_B',  ray_primaryColour(3));

    %add_feature(flowerIdx, 'colour_difference_disk',  disk_colourDiff);  %overall colour difference between UPOVs within ray florets
    add_feature(flowerIdx, 'colour_difference_trans', trans_colourDiff); 
    add_feature(flowerIdx, 'colour_difference_ray',   ray_colourDiff); 

    %add_feature(flowerIdx, 'L_trans_mean', nanmean(gerb.LAB_masked.trans(:,:,1),[1 2]));
    %add_feature(flowerIdx, 'A_trans_mean', nanmean(gerb.LAB_masked.trans(:,:,2),[1 2]));
    %add_feature(flowerIdx, 'B_trans_mean', nanmean(gerb.LAB_masked.trans(:,:,3),[1 2]));

    %add_feature(flowerIdx, 'L_ray_mean', nanmean(gerb.LAB_masked.ray(:,:,1),[1 2]));
    %add_feature(flowerIdx, 'A_ray_mean', nanmean(gerb.LAB_masked.ray(:,:,2),[1 2]));
    %add_feature(flowerIdx, 'B_ray_mean', nanmean(gerb.LAB_masked.ray(:,:,3),[1 2]));
    
    %add_feature(flowerIdx, 'L_ray_secondary', ray_secondaryColour(1));
    %add_feature(flowerIdx, 'A_ray_secondary', ray_secondaryColour(2));
    %add_feature(flowerIdx, 'B_ray_secondary', ray_secondaryColour(3));
    
    %% Correlations
    add_feature(flowerIdx, 'colour_slope_disk_L', coeffsLAB{1,1}(1));
    add_feature(flowerIdx, 'colour_slope_disk_A', coeffsLAB{2,1}(1));
    add_feature(flowerIdx, 'colour_slope_disk_B', coeffsLAB{3,1}(1));
    
    add_feature(flowerIdx, 'colour_slope_trans_L', coeffsLAB{1,2}(1));
    add_feature(flowerIdx, 'colour_slope_trans_A', coeffsLAB{2,2}(1));
    add_feature(flowerIdx, 'colour_slope_trans_B', coeffsLAB{3,2}(1));
    
    add_feature(flowerIdx, 'colour_slope_ray_L', coeffsLAB{1,3}(1));
    add_feature(flowerIdx, 'colour_slope_ray_A', coeffsLAB{2,3}(1));
    add_feature(flowerIdx, 'colour_slope_ray_B', coeffsLAB{3,3}(1));
    
    %{
    add_feature(flowerIdx, 'L_disk_intercept', coeffsLAB{1,1}(2));
    add_feature(flowerIdx, 'A_disk_intercept', coeffsLAB{2,1}(2));
    add_feature(flowerIdx, 'B_disk_intercept', coeffsLAB{3,1}(2));
    
    add_feature(flowerIdx, 'L_trans_intercept', coeffsLAB{1,2}(2));
    add_feature(flowerIdx, 'A_trans_intercept', coeffsLAB{2,2}(2));
    add_feature(flowerIdx, 'B_trans_intercept', coeffsLAB{3,2}(2));
    
    add_feature(flowerIdx, 'L_ray_intercept', coeffsLAB{1,3}(2));
    add_feature(flowerIdx, 'A_ray_intercept', coeffsLAB{2,3}(2));
    add_feature(flowerIdx, 'B_ray_intercept', coeffsLAB{3,3}(2));
    %}
    
    %% Gradients 
    %add_feature(flowerIdx, 'ray_segment_location', (ray_segment-seg(2))/(seg(3)-seg(2)));  %location of the segmentation within the ray florets
    %add_feature(flowerIdx, 'trans_segment_location', (trans_segment-seg(1))/(seg(2)-seg(1)));
    %add_feature(flowerIdx, 'disk_segment_location', disk_segment/seg(1));
    
    add_feature(flowerIdx, 'gradient_boolean_ray',   double(ray_grad_bool)); %boolean value if flower contains gradient on ray florets
    add_feature(flowerIdx, 'gradient_boolean_trans', double(trans_grad_bool)); 
    %add_feature(flowerIdx, 'gradient_boolean_disk',  double(disk_grad_bool));
    
    add_feature(flowerIdx, 'gradient_position_ray',   ray_binPsychoParams(1));   %gradient position change: threshold of psychometric function changing between the primary/secondary colours
    add_feature(flowerIdx, 'gradient_position_trans', trans_binPsychoParams(1)); %gradient position change: threshold of psychometric function changing between the primary/secondary colours
    add_feature(flowerIdx, 'gradient_shape_ray',      ray_binPsychoParams(2));     %gradient shape: absolute value of slope changing between the primary/secondary colours
    add_feature(flowerIdx, 'gradient_shape_trans',    trans_binPsychoParams(2));   %gradient shape: absolute value of slope changing between the primary/secondary colours
    
    %% Sizes
    add_feature(flowerIdx, 'size_disk', seg(1)/seg(3));
    add_feature(flowerIdx, 'size_trans', (seg(2)-seg(1))/seg(3));
    add_feature(flowerIdx, 'size_ray', (seg(3)-seg(2))/seg(3));
    %add_feature(flowerIdx, 'size_golden_ratio_sumsqs', gr_ss);

    %% GLCM / Tamura features
    add_feature(flowerIdx, 'GLCM_Tamura_coarseness',     TamuraStats(1));
    add_feature(flowerIdx, 'GLCM_Tamura_contrast',       TamuraStats(6));
    add_feature(flowerIdx, 'GLCM_Tamura_directionality', TamuraStats(5));
    
    %offset_size = 1;
    %add_feature(flowerIdx, 'GLCM_Matlab_contrast',    glcm_stats(1,offset_size));
    %add_feature(flowerIdx, 'GLCM_Matlab_correlation', glcm_stats(2,offset_size));
    %add_feature(flowerIdx, 'GLCM_Matlab_energy',      glcm_stats(3,offset_size));
    %add_feature(flowerIdx, 'GLCM_Matlab_homogeneity', glcm_stats(4,offset_size));

    %% NGTDM features
    add_feature(flowerIdx, 'NGTDM_coarseness', ngtdm_stats.coarseness);
    add_feature(flowerIdx, 'NGTDM_contrast',   ngtdm_stats.contrast);
    add_feature(flowerIdx, 'NGTDM_complexity', ngtdm_stats.complexity);
    
    %add_feature(flowerIdx, 'NGTDM_busyness', ngtdm_stats.busyness); %calculation seems off
    %add_feature(flowerIdx, 'NGTDM_strength', ngtdm_stats.strength); %poorly defined by original article
    
    %% NGTDM other version
    %add_feature(flowerIdx, 'NGTDM_other_coarseness', other_ngtdm_stats.coarseness); 
    %add_feature(flowerIdx, 'NGTDM_other_contrast', other_ngtdm_stats.contrast);
    %add_feature(flowerIdx, 'NGTDM_other_busyness', other_ngtdm_stats.busyness);
    %add_feature(flowerIdx, 'NGTDM_other_complexity', other_ngtdm_stats.complexity);
    %add_feature(flowerIdx, 'NGTDM_other_strength', other_ngtdm_stats.strength);
    
    %% Properties from flowerpoll survey
    %attribute_idxs = [1 2 3 4 8 9]; %bullseye, busyness, complexity, depth, pointiness, symmetry (appeal and interest are removed manually)
    %for i=attribute_idxs %2:length(poll.headers) %col 1 is appeal
    for i=1:length(poll.headers) 
        add_feature(flowerIdx, ['poll_' poll.headers{i}], poll_attributes(i));
    end

    fprintf('\n'); %next line for user output
end

%save output for RSA
save([dataDir 'properties_table_' imgAbbr],'imgProperties','header','labels');

end

%{
function [diffVal] = relativeIntensityJudgement(newVal, oldVal, intensityDiff)

if newVal > (oldVal + intensityDiff)
    diffVal = 1;
elseif newVal < (oldVal - intensityDiff)
    diffVal = 0;
else
    diffVal = oldVal;
end

end
%}

function add_feature(flowerIdx, label, value)
%adds feature to next column in list of features and corresponding label

global feature_count;
global imgProperties;
global header;

feature_count = feature_count+1; %increment column

imgProperties(flowerIdx, feature_count) = value;

if flowerIdx==1
    header{feature_count} = label; 
end 

end