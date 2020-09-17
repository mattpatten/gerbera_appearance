function Zr = findRadialAvgs(img,npix,xo,yo,method,omitnan)

% Adapted version of the RADIALAVG function available on Matlab Central.
%
% RADIALAVG	 Radially averaqe 2D square matrix z into m bins
%http://au.mathworks.com/matlabcentral/fileexchange/46468-radialavg-zip?focused=6478826&tab=function

% [Zr] = findRadialAvgs(img,npix,xo,yo,method,omitnan)
%
% [Zr] = findRadialAvgs(img,npix,xo,yo,method,omitnan) computes the average along the 
% radius of a unit circle inscribed in the square matrix img. The average is computed 
% in bins of size 'npix'. The radial average is not computed beyond the unit circle, 
% in the corners of the matrix img. The radial average is returned in Zr. Not a Number
% (NaN) values are excluded from the calculation. If offset values xo,yo are used, the 
% origin (0,0) of the unit circle about which the findRadialAvgs is computed is offset 
% by xo and yo relative to the origin of the unit square of the input img matrix.
%
% Example
%	N=101;
%	[X,Y] = meshgrid(-1:2/(N-1):1);
%	xo = +0.25;
%	yo = -0.25;
%	X = X-xo;
%	Y = Y-yo;
%	z = 1-sqrt(X.^2 + Y.^2);
%	m=(N-1)/2+1;
%	[Zr] = findRadialAvgs(z,m,xo,yo,'mean',1);
%	figure;plot(Zr,'.-');
%
% INPUT
% img     = square input matrix to be radially averaged
% npix    = number of pixels wide for each annulus     (MattP)
% xo      = offset of x-origin relative to unit square (DEF: 0)
% yo      = offset of y-origin relative to unit square (DEF: 0)
% method  = Calculation method. Use 'mean' or 'median' (MattP)
% omitnan = logical for removing nans during analysis  (MattP)
%
% OUTPUT
% Zr  = radial average of length m
%
% (c) 2014 David J. Fischer | fischer@shoutingman.com
% 4/4/14 DJF first working version
% 5/2/14 DJF documentation & radialavg_tester.m to demonstrate use
% radial distances r over grid of z
% 6/20/16 DJF Excludes NaN values
% 6/21/16 DJF Added origin offset
% 2/11/2018 MLP Hacked to define by number of pixels, not number of bins. 
%           Added mean/median option and handling of NaN values.

[xImg, yImg] = size(img);
%if no centre defined, assume no offset
if ~exist('xo','var')
	xo = xImg/2;
end
if ~exist('yo','var')
	yo = yImg/2;
end

N = min(size(img));
[X,Y] = meshgrid(linspace(-N/2,N/2,N)); %create grid
X = X-(xo-N/2); %define centre
Y = Y-(yo-N/2);

%this generates the number of pixels away from the centre
r = round(sqrt(X.^2+Y.^2)); %get polar distance values for each space in the matrix
ctr = r==1;

% equi-spaced points along radius which bound the bins to averaging radial values
rbins = 1:npix:N/2;
nbins = length(rbins); %last element

Zr = zeros(1,length(rbins-1)); % vector for radial average
nans = ~isnan(img); % identify NaNs in input data

% loop over the bins, except the final bin
for j=1:nbins-1
	% find all matrix locations whose radial distance is in the jth bin
	bins = r>=rbins(j) & r<rbins(j+1);
	
	% exclude data that is NaN
    if omitnan, bins = logical(bins .* nans); end
    
	% count the number of those locations
	n = sum(sum(bins));
	if n~=0
		% average the values at those binned locations
        if strcmpi(method,'mean')
            Zr(j) = sum(img(bins))/n;
        elseif strcmpi(method,'median')
            Zr(j) = median(img(bins));
        elseif strcmpi(method,'mode')
            [~, Zr(j)] = max(histcounts(img(bins),'BinWidth',1));
        elseif strcmpi(method,'sum')
            Zr(j) = sum(img(bins)) / sum(sum(bins));
        end
	else
		% special case for no bins (divide-by-zero)
		Zr(j) = NaN;
	end
end

% special case the last bin location to not average Z values for
% radial distances in the corners, beyond R=1
bins = r>=rbins(nbins) & r<=N/2;

% exclude data that is NaN
if omitnan, bins = logical(bins .* nans); end

n = sum(sum(bins));
if n~=0
	% average the values at those binned locations
    if strcmpi(method,'mean')
        Zr(nbins) = sum(img(bins))/n;
    elseif strcmpi(method,'median')
        Zr(nbins) = median(img(bins));
    elseif strcmpi(method,'mode')
        [~, Zr(nbins)] = max(histcounts(img(bins),'BinWidth',1));
    elseif strcmpi(method,'sum')
        Zr(nbins) = sum(img(bins)) / sum(sum(bins));
    end
else
	Zr(nbins) = NaN;
end
end