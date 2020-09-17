function [segmentation] = radial_segmentation(img, p, lbls)

% Performs segmentation of flower using:
%   * k-means clustering based on changes in colour, to identify the boundary between the disk and trans florets.
%   * percentiles of the kernel density estimation from the phase congruency image, to identify the boundary 
%   between the trans and ray florets (petals).
% 
% The analysis is performed radially, starting at the centre and growing in eccentricity. Due to the limited 
% number of pixels at the centre, it is recommended to remove a small number of eccentricites (~3).
%
% %%%%%%%%
% A full list of parameters and inputs are listed and described in transform_and_segment.m
% %%%%%%%%
%
% Output is a 2-column matrix of positive integers saved to a file within the image directory. 
% The first column is the radial distance from the centre to the end of the disk florets for each of the flowers 
% and the second column is the distance from the centre to the end of the trans florets.
% 
% Created by Matt Patten in Nov, 2018.

%get some properties
[xSize, ySize, nColours] = size(img.analysis);
xCtr = round(xSize/2); %flower is shifted to centre of image already
yCtr = round(ySize/2); 
kmeansColours = colormap(jet(p.num_centroids)); close; %choose colours for each cluster

%initalize
rad_median = zeros(floor(min([xSize ySize])/2),1);
rad_mean   = zeros(floor(min([xSize ySize])/2),3);

%get mesh grid (matrix of the distance of each pixel to the centre)
r = convert2polar([xSize ySize], [xCtr yCtr]);
ctr = or(r==0,r==1);

%Compute radial values
rad_median(:,1) = findRadialAvgs(img.phase, p.num_pixels, xCtr, yCtr, 'median', 0); %calculate median phase value for each eccentricity
for cc=1:nColours
    rad_mean(:,cc) = findRadialAvgs(img.analysis(:,:,cc), p.num_pixels, xCtr, yCtr, 'mean', 0); %calculate mean colour for each eccentricity
end
numBins = size(rad_mean,1);

%Cluster analysis using colour information
cluster3D = perform_kmeans(rad_mean, p.num_centroids, p.remove_pixels); %perform k-means clustering
cluster3D = cluster3D(~isnan(cluster3D)); %remove any trailing NaN's from outside the mask
numBins_afterMask = length(cluster3D);
boundaries = [0; (diff(cluster3D))]; %when there's a change in cluster label
kmeans_idxs = find(boundaries~=0)'; %get the locations of these changes

%Fit kernel density estimate (after unbinning median data)
[rad_kdeFit, phase_idxs] = kde_from_binned_data(1:numBins, rad_median, p.multiplier, p.sigValue);

%final result - combining colour segmentation (inner) and phase segmentation (outer)
segmentation = [kmeans_idxs(1) phase_idxs(end)];


%% Display results
if p.dispFig
    
    %get circles associated with the boundary indexes
    kmeans_circles = zeros(size(r)); %initialize
    for s=1:length(kmeans_idxs)
        kmeans_circles = kmeans_circles + double(r==kmeans_idxs(s));
    end
    phase_circles = (r==phase_idxs(1) | r==phase_idxs(2));
    segmentation_circles  = (r==segmentation(1)  | r==segmentation(2));
    
    %for subplot
    numRows = 3;
    numCols = 5;
    imgTitles = {'L: Luminance','A: Green to Red','B: Blue to Yellow'};
    plotTitles = {'L: Intensity by ecc','A: Intensity by ecc','B: Intensity by ecc'};
    
    f = figure;
    
    %% Top Row
    
    %draw segmentation from both colour and phase information
    subplot(numRows,numCols,1);
    imshow(imoverlay(imoverlay(img.display,ctr,'r'),segmentation_circles,'c'));
    title('Segmentation (Colour [inner], Phase [outer])');
    axis image; axis off;
    
    %histogram
    subplot(numRows, numCols, 2);
    histogram(round(rgb2gray(img.display)));
    title('Image Histogram (Greys)');
    axis square tight; box off;
    
    
    %draw image with k-means segmentation boundaries
    subplot(numRows,numCols,4);
    imshow(imoverlay(imoverlay(img.display,ctr,'r'),kmeans_circles,'g'));
    title('Segmentation (Colour)');
    axis image; axis off;
    
    %2nd row - mask
    subplot(numRows, numCols, numCols+1);
    imshow(imread([lbls.imgDir 'proc' filesep 'ctr_' lbls.image '.png']));
    title('Mask + Centre');
    axis image; axis off;
    
    %% Draw image / plot of each analysed colour channel
    
    for col=1:nColours
        subplot(numRows, numCols, numCols+1+col);
        
        %draw image
        imagesc(img.analysis(:,:,col));
        colormap(gca,'gray');
        title(imgTitles{col});
        axis image; axis off;
        
        %draw eccentricity figure, marking clusters
        subplot(numRows,numCols,numCols*2+1+col);
        hold all;
        line([kmeans_idxs; kmeans_idxs],[min(rad_mean(:,col)) max(rad_mean(:,col))+0.05],'Color','g'); %vertical line boundaries
        for i=1:numBins_afterMask
            plot(i,rad_mean(i,col),'Marker','o','MarkerEdgeColor','None','MarkerFaceColor',kmeansColours(cluster3D(i),:),'MarkerSize',3); %plot data
        end
        set(gca,'XLim',[1 numBins],'XTick',linspace(0,numBins,6));
        xlabel('Eccentricity bin'); ylabel('Mean pixel intensity'); title(plotTitles{col});
        axis square tight; box off;
    end
    
    %draw 3d clustered plot
    subplot(numRows,numCols,numCols*2+1);
    for i=numBins_afterMask:-1:1 %go backwards, for presentation reasons
        plot3(rad_mean(i,1),rad_mean(i,2),rad_mean(i,3),'Marker','o','MarkerFaceColor',kmeansColours(cluster3D(i),:),'MarkerEdgeColor','k','MarkerSize',5);
        hold all;
    end
    set(gca,'Ydir','reverse');
    xlabel('L'); ylabel('A'); zlabel('B'); %title('Mean Per Ecc (3-channel)');
    axis square; grid on;
    
    
    %% Draw phase congruency (segmentation, image and kde plot)
    
    %draw image with phase congruency segmentation boundaries
    subplot(numRows,numCols,numCols*1);
    imshow(imoverlay(imoverlay(img.display,ctr,'r'),phase_circles,'m'));
    title('Segmentation (Phase)');
    axis image; axis off;
    
    subplot(numRows,numCols,numCols*2);
    %imagesc(imoverlay(imoverlay(imgs.phase,ctr,'r'),segment_estimates,'g'));
    imagesc(img.phase);
    colormap(gca,'gray');
    title('Phase Congruency Img');
    axis image; axis off;
    
    subplot(numRows,numCols,numCols*3); hold all;
    line([phase_idxs; phase_idxs],[min(rad_median) max(rad_median)+0.05],'Color','m'); %vertical line boundaries
    plot(1:numBins_afterMask,rad_kdeFit(1:numBins_afterMask),'LineStyle','-','Marker','None','Color','k','LineWidth',1.5); %plot kde
    for i=1:numBins_afterMask
        plot(i, rad_median(i),'Marker','o','MarkerEdgeColor','None','MarkerFaceColor',kmeansColours(cluster3D(i),:),'MarkerSize',3); %plot data
    end
    
    %axis properties
    set(gca,'XLim',[1 numBins_afterMask],'XTick',linspace(0,numBins,6));
    xlabel('Eccentricity bin'); ylabel('Median pixel intensity'); title('G: Intensity of Phase Cong by ecc');
    axis square tight; box off;
    
    %Set position on screen and save
    set(gcf, 'Position', [1, 1, 1600, 1200]); %set size on screen
    saveas(f,[lbls.dir 'Segmentation_' lbls.image '_' lbls.analysis '.png']);
    close;
end

end