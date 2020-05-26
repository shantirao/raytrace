function psf = pupilPSF(opd, amplitude, bands, intensities, focalLength, npix, pupilSamplingDistance, pixelSize)
% this method was developed by Andy Kee. The built-in FFT function
% uniformly samples frequencies with interval 1/n. But we don't want that!
% We only want the intervals that correspond to pixels. Using FFT to
% calculate a PSF means padding the pupil to a big matrix, doing too many
% computations, and then downsampling to throw away most of them. Instead,
% Andy's DFT2 function can do the job much faster
rays = [];
if iscell(opd)
    opd = opd{end};
end
if isstruct(opd)
    rays = opd;
    [ opd, amplitude] = pupilWFE(opd);
elseif nargin < 2
    amplitude = ~isnan(opd);
end
if nargin < 3 
    bands = [];
end
if isempty(bands)
    bands = .0006328;
end
if nargin < 4
    intensities = [];
end
if isempty(intensities)
   intensities = ones(numel(bands),1)/numel(bands);
end
if nargin < 5
    npix = 100;
end
if nargin < 6
     if isstruct(rays) && isfield(rays.samplingDistance)
        samplingDistance = rays.samplingDistance;
     else
         pupilSamplingDistance = 1; % mm
     end
end
if nargin < 7
    if isstruct(rays) && isfield(rays.pixelSize)
        pixelSize = rays.pixelSize;
    else
        pixelSize = 0.010; % 10 micron
    end
end
if nargin < 8
    focalLength = 100; %makes for a nice picture
end

sampling = pupilSamplingDistance .* pixelSize / focalLength;

M = npix(1);
N = npix(end);
[m,n] = size(opd);

X = (-(n-1)/2:(n-1)/2)';
Y = (-(m-1)/2:(m-1)/2)';
U = (-(N-1)/2:(N-1)/2)';
V = (-(M-1)/2:(M-1)/2)';
% X = (-(n)/2:(n)/2-1)';
% Y = (-(m)/2:(m)/2-1)';
% U = (-(N)/2:(N)/2-1)';
% V = (-(M)/2:(M)/2-1)';

psf = zeros(M,N);

for k = 1:length(bands)
    lambda = bands(k);
    alpha = sampling./lambda;
    psf = psf + intensities(k) * abs(exp(-2*pi*1j*alpha(end)*V*Y') * (amplitude .* exp(1j*2*pi/lambda*opd)) * exp(-2*pi*1j*alpha(1)*X*U')).^2;
end

end