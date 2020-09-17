function show_image_subset(imgAbbr, labels)

% send in a bunch of labels (filenames) and this will display the bunch for you
%
% Created by Matt Patten
% Created in February 2020


imgDir = get_dir(imgAbbr, 'img');

cols = 10;
rows = ceil(length(labels)/cols);
figure;

for flowerIdx=1:length(labels)
    gerb = imread([imgDir  labels{flowerIdx}  '.png']); %load image
    subplot(rows,cols,flowerIdx);
    imshow(gerb);
    title(labels{flowerIdx});
    axis image;
end

set(gcf, 'Position', [1, 1, 1600, 1000]); %set size on screen

end
