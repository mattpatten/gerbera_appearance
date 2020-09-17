function fixSegmentationIncremental(imgAbbr)

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

imageIdx = [];
while isempty(imageIdx)
    imageNumStart = inputPosIntQ('What flower would you like to start manual correction at (or press q to quit)?');
    if imageNumStart==0, return; end %exit function if user has pressed q (aka 0)
    imageIdx = find(strcmpi(labels,[imgAbbr sprintf('%04i',imageNumStart)]));
    if isempty(imageIdx), disp('Flower index does not exist in this dataset. Please check your input.'); end
end

while imageIdx<=length(labels)
    
    disp(['Loading flower ' labels{imageIdx} '...'])
    
    %load image for this flower
    gerb.RGB  = imread([imgDir  labels{imageIdx}  '.png']); %load image
    seg = segmentation(imageIdx,:);
    
    %display automated segmentation on image
    f = drawFigure(gerb.RGB, r, seg(1), seg(2));
    
    %check if there are more flowers to edit
    proceed = inputYesNoQ('Is this a satisfactory segmentation for this flower (y+/n- or q)? ');
    close(f);
    if strcmpi('y',proceed) || strcmpi(proceed,'+')
        imageIdx = imageIdx+1;
        
    elseif strcmpi('q',proceed)
        break;
        
    else
        %manually choose segmentation region
        for i=1:2 %for each circle
            
            %display automated segmentation on image
            f = drawFigure(gerb.RGB, r, seg(i));
            
            userHappy = false;
            firstTime = true;
            while ~userHappy
                
                if ~firstTime
                    %update segmentation
                    rad = inputPosIntQ(instructions{i});
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
                resp = inputYesNoQ('Is this satisfactory (y+/n-)? ');
                if strcmpi(resp,'y') || strcmpi(resp,'+')
                    userHappy = true;
                end
            end
            close(f);
        end
        
        segmentation(imageIdx,:) = seg;
        
        %show completed segmentation
        f = drawFigure(gerb.RGB, r, seg(1), seg(2));
        
        %check if there are more flowers to edit
        proceed = inputYesNoQ('Would you like to process another flower (y+/n-)? ');
        if strcmpi('n',proceed) || strcmpi(proceed,'-') || strcmpi(proceed,'q')
            close(f);
            break;
        end

        imageIdx = imageIdx+1; %iterate
        close(f);
        
    end
end

%Careful - this will overwrite existing files
save([imgDir 'segmentation_manual.mat'],'segmentation');

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


function [val] = inputPosIntQ(inString)
%get input and ensure is positive integer

suitableInput = false;

while ~suitableInput
    val = input(inString,'s');
    if strcmpi('q',val)
        val = 0;
        suitableInput = true;
    else
        val = str2double(val);
        if rem(val,1)==0 && val>0
            suitableInput = true;
        else
            disp('Invalid input: Positive integers only.');
        end
    end
end

end


function [val] = inputYesNoQ(inString)
%get input and ensure value is y or n

suitableInput = false;

while ~suitableInput
    val = input(inString,'s');
    if strcmpi('n',val) || strcmpi('y',val) || strcmpi('+',val) || strcmpi('-',val) || strcmpi('q',val)
        suitableInput = true;
    else
        disp('Invalid input. y/+/n/- or q (save and quit) only.');
    end
end

end