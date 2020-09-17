function [thr] = decide_rg_threshold(img)

%Inputs greyscale image and decides threshold for region growing algorithm based on
%mean intensity values (i.e., images with mean intensity that is closer to white 
%should have smaller threshold than those that are quite distinct)
%
%This is important because the algorithm will remove white flowers if the threshold is
%too high, and it will do a very botch job otherwise, with lots of squares and jagged 
%edges around the outside of the petals if the threshold is too low. So the idea is to 
%keep the threshold high unless this would cause a problem (i.e., white flowers).
%
%Created by Matt Patten in Dec 2018

standard_threshold = 20;
img_mean = mean(mean(img));

%thr = 4 * log(255 - img_mean); %increases too sharply, meaning mask for white flowers leaks into petals
%linearly increases away from white and then caps at specified value. NB: For -(A/B), at B points away from white, the threshold will be A
thr = min([standard_threshold, -(5/15).*(img_mean-255)]);  %-5/15 <-- standard for pref/united/pin data

disp([',  Image mean from white: ' num2str(255 - img_mean) ' Threshold: ' num2str(thr)]);

end