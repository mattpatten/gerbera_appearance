function run_pipeline

% Everything that needs to be done from once an image is downloaded to the RSA analysis.
%
% More specifically:
%   (i)    Rename the image to a standardized format (img01.png).
%   (ii)   Finding the flower face in the image (via region growing algorithm from the image corners), and generating a mask for 
%          regions outside of this area.
%   (iii)  Identify flower centre (fminsearch of circular area within mask), pad/crop surround to make image a perfect square 
%          and rescale image to specific dimensions (e.g., 300x300).
%   (iv)   Perform transformations of images as required (L*a*b*/grey colourspaces, applying mask, finding pixels along mask border)
%   (v)    Perform Kovesi's phase congruency algorithm to identify stimulus edges and high frequency areas.
%   (vi)   Perform k-means clustering on the L*a*b* colour channels to identify eccentricities with rapid colour changes.
%   (vii)  Segment flower based on the above two algorithms, saving a mat file with the eccentricities at the edge of the disk and
%          trans florets for each flower.
%   (viii) Allow manual corrections to segmented flowers and save images that display the segmentation values.
%   (ix)   Generates values for ~45 properties of the flower, storing it in a table and displaying images in the dataset that are 
%          the highest and lowest on each of these properties, and combining results for all flowers into histograms.
%   (x)    Uses the values across all of these properties to group the flower set into a specific number of clusters.
%   (xi)   Performs a regression analysis on these properties in order to explain behavioural ratings measured elsewhere (see flwrpoll)
%          and identify the attributes that contribute to human preference.
%
% Created by Matt Patten on 27 Nov 2018

%parameters
imgAbbr = 'ABC'; %Abbreviation used for file names and structure, indicating the dataset to be used. So far: 'pref','pres','orig','pin','united','fc','buy'
pixSize = 300;
sourceDir = ''; % relative or hard-coded path to stimulus image folder - required only for original non-processed images

%add libraries/functions
addpath(genpath([fileparts(mfilename('fullpath')) filesep 'functions' filesep]));      %look for sub-directory from current m-file called 'functions' and add any sub-directories
addpath(genpath([fileparts(mfilename('fullpath')) filesep 'segmentation' filesep]));   %look for sub-directory from current m-file called 'segmentation' and add any sub-directories
addpath(genpath([fileparts(mfilename('fullpath')) filesep 'img_properties' filesep])); %look for sub-directory from current m-file called 'img_properties'. Don't add any of its sub-directories
addpath([fileparts(mfilename('fullpath')) filesep 'flowerpoll' filesep]);              %look for sub-directory from current m-file called 'flowerpoll'. Don't add any of its sub-directories
addpath([fileparts(mfilename('fullpath')) filesep 'regression' filesep]);              %look for sub-directory from current m-file called 'regression'. Don't add any of its sub-directories

%Pipeline processes
rename_images(imgAbbr, sourceDir); %load all images in specified directory and rename them to standard name with incremental numbers - be wary using this as don't want to rename flowerIDs shown on the website
create_masks(imgAbbr,pixSize);       %generate masks for these images (if not already provided)
remove_bad_images(imgAbbr);          %remove any images where masking failed (messy background, too similar to background)
resize_images(imgAbbr,pixSize);      %find centre and resize images to a specific, standardized size.
transform_and_segment(imgAbbr);      %NB: some parameters are located within this
%fixSegmentationIncremental(imgAbbr); %Manually apply corrections to any segmentations that did not go as intended, going through each flower one-by-one
fixSegmentationManually(imgAbbr);   %Manually apply corrections to any segmentations that did not go as intended, choosing specific flowers to fix
saveSegmentationComparison(imgAbbr); %save all (both auto and those manually corrected) segmentation images in their own folder
saveSegmentationFinal(imgAbbr);      %save the final segmentation images in their own folder - QA / sanity check - only takes a second to run
get_flower_properties(imgAbbr);      %Extracts various statistical properties of the flower into a table for analysis
generate_examplar_images(imgAbbr);   %loads up the property table and looks at the lowest and highest on each score
generate_histograms(imgAbbr);
cluster_imageset(imgAbbr, 10, 90);   %runs dimensionality reduction and pca on the properties table to divide the image set into x groups when explaining y% of the variance.
regressionAnalysis('appeal',1);

end