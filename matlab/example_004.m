% Example 004: chromatic ray reflects from diffraction grating

% [start, finish] = importVV('example_001_VV.txt');

%% setup
% prism front surface
s1=struct();
s1.position = [0,0,60];
s1.direction=normr([0,0,-1]);
s1.type = 'reflect';
s1.grating = cross([1,0,0],s1.direction)/1e-2; %ruling direction divided by spacing, 10 microns
s1.order = -1;

% Stop
s3=struct();
s3.position = [0,0,0];
s3.direction=[0,0,1];

%% trace 1e3 times in 2 seconds
N = 12;
source = struct();
source.position=repmat([0,0,0],N,1);
source.direction=repmat(normr([0,0,1]),N,1);
source.valid=true(N,1);
source.wavelength=.000400+.000025*(1:N)';
source.units='mm';

surfaces = {s1, s3};

tic
trace = raytrace(source, surfaces);
toc
rays=trace{end};
%%
clf;
viewaxes = [0 0 1;0 1 0]';
plotSideView(trace,viewaxes,'LineWidth',2);
colororder(gca,wavelengthToRGB(source.wavelength));
hold on;
plotLine(s1.position*viewaxes,(s1.position + 10*cross(s1.direction,[1,0,0]))*viewaxes,'k');
axis image; xlabel('z'); ylabel('y');
hold off;

%% compare

y = rays.position(:,2);
hypot = rays.opl - 60;
sinRefracted = y./hypot;
lambdaOverD = rays.wavelength / 1e-2; % wavelength / grating pitch
disp(sprintf('Verify: %d should be close to 0',max(abs(sinRefracted - lambdaOverD))));
