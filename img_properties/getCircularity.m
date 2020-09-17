function getCircularity(img, imgAbbr, label)

% This computes 1) the phase congruent polar image and takes the outside edge
%               2) the edge of the mask of the polar-transformed image 
% and estimates the derivative of each point, displaying this in a histograph.
%
% Inputs:
%    img   - A struct for a series of images.
%    imgAbbr - 
%    labels - Strings providing a name/label for each flower.
%
% Output:
%     A series of graphs.
%
% Created by Matt Patten in March, 2019.


% Initialize
outputDir = get_dir(imgAbbr,'output');
circDir = [outputDir 'polarCircularity' filesep]; 
if ~exist(circDir,'dir'), mkdir(circDir); end
[xSize, ySize, ~] = size(img.analysis); %get image properties


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Phase congruency     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define parameters
p.pc.k = 20; %checked between 5-20, 10-20 is the sweet spot
p.pc.nscale = 2; %checked between 2-6: 2 is the best, 3 isn't so bad either (int only)
p.pc.minWavelength = 4; %best to leave centre but get rid of rest

%perform phase congruency analysis
img.phaseCong = phasecong3(img.analysis/255, 'k', p.pc.k, 'nscale',p.pc.nscale, 'minWaveLength', p.pc.minWavelength);
img.phaseBinary = imbinarize(img.phaseCong,'adaptive');
img.phaseBinary([1 end],:) = 0; %edges have phase shift, but is not relevant in our case


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calc angularity      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgeImg = zeros(xSize,ySize);
[~,j] = max(flipud(img.phaseBinary));
edges_pc = xSize-j+1;
edgeImg(sub2ind([xSize ySize],edges_pc,1:ySize)) = 1;
img.pcEdges = edgeImg;


%generate map of edges of flower after polar transform
img.maskEdge = [zeros(1,ySize); abs(diff(img.mask))];
[edges_mask,edge_vals] = find(img.maskEdge~=0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Plotting         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%draw figure
f = figure; 
set(gcf,'visible','off'); %stops window popping up every split second so we can go do other things
rows = 2; cols = 5;

%phase cong
subplot(rows,cols,1); imshow(img.regularDisplay);   title('Original');
subplot(rows,cols,2); imshow(img.phaseCong);        title('Phase Congruency');
subplot(rows,cols,3); imshow(img.phaseBinary);      title('Binarized PC');
subplot(rows,cols,4); imshow(img.pcEdges);          title('Phase Cong Edge');
data = []; data = abs(diff(edges_pc));
subplot(rows,cols,5); hist(data,0:max(data)); title(['PC Hist (' num2str(sum(data)) ')']); xlim([-0.5 11]); xticks(0:10); box off;

subplot(rows,cols,6); imshow(img.polarDisplay/255); title('Polar Image');
subplot(rows,cols,8); imshow(img.mask);             title('Mask');
subplot(rows,cols,9); imshow(img.maskEdge);         title('Mask Edge');
data = []; data = abs(diff(edges_mask));
subplot(rows,cols,10); hist(data,0:max(data)); title(['Mask Hist (' num2str(sum(data)) ')']); xlim([-0.5 11]); xticks(0:10); box off;

set(gcf, 'Position', [1, 1, cols*300, rows*400]); %set size on screen
saveas(f,[circDir 'circularity_' label '.png']);
close;

end