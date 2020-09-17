function saveSegmentationComparison(imgAbbr)

% This function writes each images twice, with the auto and manual segmentations side by side,
% with the overlay of any changes presented in a different colour. It also outputs how many manual 
% changes have taken place.
%
% Inputs:
%   imgAbbr - The abbreviation of the dataset we're working with
%
% Output:
%   A directory filled with sets of images with auto and manual segmentations.
%
% Created by Matt Patten in January, 2019.


[imgDir, outputDir] = get_dir(imgAbbr,'img','output');
segDir = [outputDir 'segmentation_comparison' filesep];
if ~exist(segDir,'dir'), mkdir(segDir); end %if directory doesn't exist, create it

%load labels (filenames) and segmentation details
load([imgDir 'labels_' imgAbbr '.mat'],'labels');
load([imgDir 'segmentation_manual.mat'],'segmentation'); %load manual segmentations (if file is there)
segmentation_manual = segmentation; clear segmentation;

load([imgDir 'segmentation.mat'],'segmentation');     %otherwise load automatically-generated version
segmentation_auto = segmentation;

%count differences
num_edits = sum(sum((segmentation_auto - segmentation_manual)~=0));
num_segs = numel(segmentation_auto);
percent_edits = num_edits/num_segs*100;

disp(['There have been ' num2str(num_edits) ' edits of ' num2str(num_segs) ' total segmentations (' num2str(percent_edits) '%)']);

for flowerIdx = 1:length(labels)

    if mod(flowerIdx,100)==0
        disp(['Processed ' num2str(flowerIdx) ' flowers...']);
    end
    
    %load image for this flower
    filename = labels{flowerIdx};
    outputFile = [segDir filename '.png'];
    gerb.RGB  = imread([imgDir filename '.png']); %load image
    seg_man  = segmentation_manual(flowerIdx,:);
    seg_auto = segmentation_auto(flowerIdx,:);
    
    if seg_man(1)~=seg_auto(1)
        inner_colour='m';
        inner_text = 'red';
    else
        inner_colour='c';
        inner_text = 'black';
    end
    
    if seg_man(2)~=seg_auto(2)
        outer_colour='m';
        outer_text = 'red';
    else
        outer_colour='c';
        outer_text = 'black';
    end
    
    %get properties
    [xSize, ySize, ~] = size(gerb.RGB);
    r = convert2polar([xSize ySize]);

    circles_auto = or(r==seg_auto(1), r==(seg_auto(2)));
    circle_man_inner = r==seg_man(1);
    circle_man_outer = r==seg_man(2);
    
    f = figure; 
    set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
    subplot(1,2,1);
    imshow(imoverlay(gerb.RGB,circles_auto,'c'));
    title(['Auto: ' num2str(seg_auto(1)) ' & ' num2str(seg_auto(2))]);
    
    subplot(1,2,2);
    imshow(imoverlay(imoverlay(gerb.RGB,circle_man_inner,inner_colour),circle_man_outer,outer_colour));
    title(['Manual: {\color{' inner_text '}' num2str(seg_man(1)) '} & \color{' outer_text '}' num2str(seg_man(2))]);
    
    %imwrite(imoverlay(gerb.RGB,circles,'c'),outputFile,'png');

    %some tweaking of overall image properties
    set(gcf, 'Position', [1, 1, 900, 600]); %set size on screen
    saveas(f,outputFile);
    close;
end
end