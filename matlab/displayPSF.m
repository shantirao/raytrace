% -------------------------------------------------------------------------------
% [psf] = DisplayPSF(opd, ValidIdx, units, bandpass);
% -------------------------------------------------------------------------------
%
% Calls MakeBroadbandPSF and displays result, with assumptions appropriate
% for the OSML example
%
% INPUTS:
%
%      opd : pupil optical path difference function 
%           (same units as wavelenth and pix_size)
%
%      ValidIdx : pupil shape
%
%      units : 
% bandpass - spectral range and relative response of detection process
%
%            (e.g.)       Wavelength           Intensity Xmission
%                   ----------------------------------------------
%                   bandpass (1,1) = 631e-9    bandpass(1,2) = 0.5
%                   bandpass (2,1) = 633e-9    bandpass(2,2) = 1.0
%                   bandpass (3,1) = 635e-9    bandpass(3,2) = 0.5
%

function [psf, scale, bandpass] = displayPSF(opd, mask, units, bandpass, FNum,pixelSize)
%DisplayPSF Helper function for plotting a point-spread function

% addpath('../../psf/image-psf-generation');
% addpath('../../psf/image-generation');
% addpath('../../psf/util');
% addpath('../../psf/zernike');

if nargin < 3
    units = 1e-3;
end
if nargin < 4
    bandpass = [632.8e-9 / units, 1]; %units are mm
end
if nargin < 5
    FNum = 1;
end
if nargin<6
    pixelSize = 1e-6 / units;
end
s = size(mask);
pupil = blkdiag(zeros(s),mask,zeros(s));
opd = blkdiag(zeros(s),opd,zeros(s));

Nyquist             = 2; %max(2,2^nextpow2(1024/length(opd)));
dfs_tilt            = 0;
%Assume a pixel size p, and (2 lambda/D * focalLength = p)
%or 2 lambda fNum = p
%lambda is the diffraction-limited wavelength for that pixel size
lambda = max(pixelSize / 2 / FNum, min(bandpass(:,1)));
sampling            = [Nyquist lambda];
[psf,bandpass]      = MakeBroadbandPSF(pupil, opd, bandpass, sampling, dfs_tilt);        

s = length(psf); %psf is always a square

% pLength = pixel * (s - 1); % physical length
% length per PSF pixel
scale = pixelSize /Nyquist;
% scale = bandpass(1,1)*Nyquist/s; % convert length to frequency

zoom = 4; %ceil(Nyquist/4); %2; %2*Nyquist;
w = ceil(1/2+s*(1-1/zoom)/2):floor(1/2+s*(1+1/zoom)/2);
a = numel(w);
x = scale*[-a/2, a/2]; 

% lpsf = (psf(w,w));
% mval = max(max(lpsf));
% imagesc(x,x,lpsf,[mval-3,mval]);
imagesc(x,x,psf(w,w));

axis image; 
% plot(psf2(round(npix/2),:));

% rmpath('../../psf/image-psf-generation');
% rmpath('../../psf/image-generation');
% rmpath('../../psf/util');
% rmpath('../../psf/zernike');
end

