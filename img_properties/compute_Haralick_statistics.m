function [stats] = compute_Haralick_statistics(imgs, p, analysisName, dispFig)

% This measures how often a pixel with the intensity (grey-level) value i occurs in a specific spatial relationship to a pixel with the value j. It is 
% set up to be a horizontal spatial relationship (i.e., a pixel to the right), but can be adjusted to be in any direction and any number of pixels away.
%
% It uses this grey-level co-occurrence matrix to compute several statistical properties of the image, such as contrast and homogeneity, as based on the
% equations presented in Haralick, Shanmugam & Dinstein (1973), Haralick (1979), Soh & Tsatsoulis (1999) and Clausi (2002).
% 
% Algorithm parameters:
%     p.bins         - The number of grey levels to consider. e.g., 32 bins means every 8 grey values are lumped into a single bin (256/32=8 values per bin).
%     p.offsetValues - Distance to measure between pixels. 1 indicates neighbouring pixels. 4 means there are three pixels in between the two pixels of interest.
%     p.symmetry     - true/false. Whether we consider an inverse relationship. i.e., [0,1] as well as [1,0]
%     p.grayLimits   - Image instensity values to scale between (i.e., min and max for greyscale values to consider, the bins will be divided within these values). Default is all: [0 255]
%
% Inputs:
%     imgs.analysis  - A series of masked, greyscale images in the order of (x, y, 1(=colour), image_index)
%     imgs.display   - For display purposes only, the original colour image.
%     p              - Parameter struct as described above.
%     analysisName   - Label for saving output files.
%     dispFig        - Whether exemplar images for each statistic should be displayed and/or saved
%
% Ouput:
%     stats          - A struct containing 22 different statistics generated from the GLCM, each struct having one statistic per image.
%
% Created by Matt Patten
% Created in February, 2019


%get properties / initialize
nImages = size(imgs.analysis,4);
glcm = zeros(p.bins,p.bins,nImages);

%% Compute GLCM
warning off; %Matlab shoots back warning that it doesn't consider NaN entries in the GLCM..... many many times
for flowerIdx = 1:nImages
    offset = [0 p.offsetValue]; %spatial relationship to the pixel being measured: [row col] e.g., [0 1] would be to the right from the reference pixel; [1 1] diagonally up-right
    glcm(:,:,flowerIdx) = graycomatrix(squeeze(imgs.analysis(:,:,1,flowerIdx)),'NumLevels',p.bins,'Offset',offset,'Symmetric',p.symmetry,'GrayLimits',p.grayLimits);
end
warning on; %not to silence all warnings - just this one

%calculate statistics according to Haralick 1979 (and others)
stats = GLCM_Features1(glcm);

%Display and save top and bottom extremes of each category/score
if dispFig
    displayExemplarImages(stats.autoc, imgs.display, 'Haralick', [analysisName '_autocorrelation']); %if the image is shifted and compared, how much overlap there would still be (and if it's periodic)
    displayExemplarImages(stats.contr, imgs.display, 'Haralick', [analysisName '_contrast']); %the amount of local variations present in an image
    displayExemplarImages(stats.corrm, imgs.display, 'Haralick', [analysisName '_correlation_matlab']);
    displayExemplarImages(stats.corrp, imgs.display, 'Haralick', [analysisName '_correlation_haralick']); %grey-tone linear dependencies
    displayExemplarImages(stats.cprom, imgs.display, 'Haralick', [analysisName '_cluster_promimence']);
    displayExemplarImages(stats.cshad, imgs.display, 'Haralick', [analysisName '_cluster_shade']);
    displayExemplarImages(stats.dissi, imgs.display, 'Haralick', [analysisName '_dissimilarity']);
    displayExemplarImages(stats.energ, imgs.display, 'Haralick', [analysisName '_energy_matlab']); % a.k.a. angular second-moment
    displayExemplarImages(stats.entro, imgs.display, 'Haralick', [analysisName '_entropy']); % "one might expect ... higher values for more complex images"
    displayExemplarImages(stats.homom, imgs.display, 'Haralick', [analysisName '_homogeneity_matlab']);
    displayExemplarImages(stats.homop, imgs.display, 'Haralick', [analysisName '_homogeneity_haralick']); %angular second-moment feature: in a homogeneous image, there ar every few dominant grey-tone transitions
    displayExemplarImages(stats.maxpr, imgs.display, 'Haralick', [analysisName '_maximum_probability']);
    displayExemplarImages(stats.sosvh, imgs.display, 'Haralick', [analysisName '_sum_of_squares_variance']);
    displayExemplarImages(stats.savgh, imgs.display, 'Haralick', [analysisName '_sum_average']);
    displayExemplarImages(stats.svarh, imgs.display, 'Haralick', [analysisName '_sum_variance']);
    displayExemplarImages(stats.senth, imgs.display, 'Haralick', [analysisName '_sum_entropy']);
    displayExemplarImages(stats.dvarh, imgs.display, 'Haralick', [analysisName '_difference_variance']);
    displayExemplarImages(stats.denth, imgs.display, 'Haralick', [analysisName '_difference_entropy']);
    displayExemplarImages(stats.inf1h, imgs.display, 'Haralick', [analysisName '_information_measure_of_correlation1']);
    displayExemplarImages(stats.inf2h, imgs.display, 'Haralick', [analysisName '_information_measure_of_correlation2']);
    displayExemplarImages(stats.indnc, imgs.display, 'Haralick', [analysisName '_inverse_difference_normalized']);
    displayExemplarImages(stats.idmnc, imgs.display, 'Haralick', [analysisName '_inverse_difference_moment_normalized']);
end

%clear some duplicate or not useful entries
clear stats.corrm; %matlab's correlation computation - duplicates .corrp (Haralick's correlation)
clear stats.energ; %matlab's energy computation - (almost) reverse of .entro: Haralick's entropy.
clear stats.homom; %matlab's homogeneity computation - duplicates .homop: Haralick's homogeneity.
clear stats.savgh; %sum average - is basically overall/average luminance (already captured in our L*a*b* measure).
clear stats.svarh; %sum variance - duplicates .sosvh: sum of squares variance.

end