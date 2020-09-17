function fixSegmentationManually(imgAbbr)

% For flowers that did not segment correctly using the automated procedure, this function
% will allow you to manually correct the location of each segmentation.
%
% NB: Only load this after automatic segmentation has been completed.
%
% Input:
%   imgAbbr - the abbreviation of the database of images we are analysing
%
% Output:
%   An edited version of the mat file storing the location of the segmentation
%
% Created by Matt Patten in January, 2019

warning off;

[imgDir] = get_dir(imgAbbr,'img');

%load image and segmentation details
load([imgDir 'labels_' imgAbbr '.mat'],'labels');

try 
    load([imgDir 'segmentation_manual.mat'],'segmentation'); %load manual segmentations (if file is there)
    disp('Loading most recent segmentation file with some manual corrections.');
catch
    load([imgDir 'segmentation.mat'],'segmentation');     %otherwise load automatically-generated version
    disp('Loading original (auto) segmentation file.');
end

%load first image as a dummy to extract relevant properties
tmp = imread([imgDir  labels{1}  '.png']);

[xSize, ySize, ~, ~] = size(tmp);
r = convert2polar([xSize ySize]);

instructions{1} = ['Please select radius (in pixels) of the INNER segmentation (between disk and trans florets: 1-' num2str(min([xSize ySize])) '): '];
instructions{2} = ['Please select radius (in pixels) of the OUTER segmentation (between trans and ray florets: 1-'  num2str(min([xSize ySize])) '): '];

clear tmp;

%check if user wants to manually fix segmentations
proceed = inputYesNo('Would you like to manually fix segmentation on any of the flowers (y+/n-)? ');
if strcmpi('n',proceed) || strcmpi(proceed,'-')
    moreFlowersToProcess = false;
else
    moreFlowersToProcess = true;
end

while moreFlowersToProcess

    success=0;
    while ~success
        %choose flower to correct
        imageNum = inputPosInt('Which flower number would you like to correct? '); %get input and ensure value is positive integer
        idx = find(strcmpi(labels,[imgAbbr sprintf('%04i',imageNum)])); %get position of image in matrix
        seg = segmentation(idx,:);

        %load image and mask for this flower
        try 
            gerb.RGB  = imread([imgDir  labels{idx}  '.png']); %load image
            success = 1;
        catch
            disp('Flower index does not exist in this dataset. Please check your input.');
            success = 0;
        end
    end
    
    %manually choose segmentation region
    for i=1:2 %for each circle
        
        %display automated segmentation on image
        f = drawFigure(gerb.RGB, r, seg(i));
        
        userHappy = false;
        firstTime = true;
        while ~userHappy

            if ~firstTime
                %update segmentation
                rad = inputPosInt(instructions{i});
                if rad < min([xSize ySize]) %save only if within acceptable range
                    seg(i) = rad;
                else
                    disp('Error: Radius too large for image dimensions.');
                end
                
                %display current segmentation on image
                close(f);
                f = drawFigure(gerb.RGB, r, seg(i));
            end
            firstTime=false;
            
            %check satisfaction of new boundary
            resp = inputYesNo('Is this satisfactory (y+/n-)? ');
            if strcmpi(resp,'y') || strcmpi(resp,'+')
                userHappy = true;
            end
        end
        close(f);
    end
    
    segmentation(idx,:) = seg;
    
    %show completed segmentation
    f = drawFigure(gerb.RGB, r, seg(1), seg(2));
    
    %check if there are more flowers to edit
    proceed = inputYesNo('Would you like to process another flower (y+/n-)? ');
    if strcmpi('n',proceed) || strcmpi(proceed,'-')
        moreFlowersToProcess = false;
    end
    
    close(f);

    %Careful - this will overwrite existing files
    save([imgDir 'segmentation_manual.mat'],'segmentation');
end

end



function [f] = drawFigure(img, r, seg1, seg2)
%draw image to window

if nargin<4, seg2=[]; end
if nargin<3, error('Not enough arguments'); end

f = figure;
set(f,'WindowStyle','docked') %dock the figure to the command window
if isempty(seg2)
    imshow(imoverlay(img,r==seg1,'c'));
    title(['Current Location: ' num2str(seg1)]);
else
    imshow(imoverlay(img,or(r==seg1,r==seg2),'c'));
    title(['Locations: ' num2str(seg1) ' & ' num2str(seg2)]);
end
commandwindow; %moves cursor back to command window instead of keeping at figure

end


function [val] = inputPosInt(inString)
%get input and ensure is positive integer

suitableInput = false;

while ~suitableInput
    val = str2double(input(inString,'s'));
    if rem(val,1)==0 && val>0
        suitableInput = true;
    else
        disp('Invalid input: Positive integers only.');
    end
end

end


function [val] = inputYesNo(inString)
%get input and ensure value is y or n

suitableInput = false;

while ~suitableInput
    val = input(inString,'s');
    if strcmpi('n',val) || strcmpi('y',val) || strcmpi('+',val) || strcmpi('-',val)
        suitableInput = true;
    else
        disp('Invalid input. y/+/n/- only.');
    end
end

end