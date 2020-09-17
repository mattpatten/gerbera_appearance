function [flag, percent_same] = check_for_natural_image(img, cols, threshold)

% Checks a set of the border/boundary pixels around the edge of the image and how much variation in 
% greyscale values there is. If there is a lot of variance, it is likely to be an image with a natural 
% background (rather than stimulus-controlled plain white backdrop, etc) and returns true.
% e.g., At least 20% of pixels in the background have to be of (close to) the same colour.
%
% Inputs:
%    img          - A 2D image matrix.
%    cols         - The number of columns around the border to examine, i.e., width, in pixels.
%    threshold    - The maximum amount of similarity (as a percentage, 0-100) that the border region 
%                   can contain for it to be considered a natural image (or minimum amount for it to be
%                   a non-natural background).
%
% Output:
%    flag         - A flag to indicate whether the image has a natural background (true) or is a monotone
%                   (false).
%    percent_same - The amount of similiarty in the border strip as a percentage of all pixels in that region.
%
% Created by Matt Patten in Dec 2018


numBins = 64; %break up grey values 0-255 into this many bins

%get border strips
edges = [img(:,1:cols);                                 ... %left columns
         img(:,(end-cols+1):end);                       ... %right columns
         img(1:cols,          (cols+1):(end-cols-1))';  ... %top row
         img((end-cols+1):end,(cols+1):(end-cols-1))'];     %bottom row
     
%flatten into single column
flat_edges = reshape(edges,numel(edges),1);

%count data, but not all 0-255 values - limit to a certain number of bins.
freq = histcounts(flat_edges,numBins); 

%find proportion of maximum bin compared to all others
percent_same = max(freq) / sum(freq) * 100; 

%Set flag, output to user
if percent_same < threshold
    flag = true;
    %disp([' REMOVED! Background similarity too low: ' num2str(percent_same)]);
else
    flag = false;
    %disp([' Background similarity: ' num2str(percent_same)]);
end

end