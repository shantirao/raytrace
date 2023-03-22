% Example 110: ordinary spherical lens, sligntly misaligned

% [start, finish] = importVV('example_200_VV.txt');

%% setup
% lens front surface
s1=struct();
s1.position = [0,0,10.427852];
s1.direction=normr([0.001,0,-1]);
s1.curvature = -1/50;  % 1 / radius of curvature, negative because it points away from the direction and it's convex
s1.type = 'refract';
s1.index = 1.518522387620793;
% direction points toward the oncoming rays, and curvature is negative,
% meaning the surface we use is concave. The sphere curves around on
% itself, right? This sign convention tells the ray propagator to use the
% intersection that falls on the outside of the sphere, not the the one
% on the interior surface of the sphere. The surface direction wants to
% point toward the oncoming beam, and cuy = 1/R is negative for convex.

% lens back surface, flat
s2=struct();
s2.position = s1.position + [0,0,10];
s2.direction=normr([0,.003,-1]);
s2.type = 'refract';
s2.index = 1;

% Stop
s3=struct();
s3.position = s1.position + [0,0,99];
s3.direction=[0,0,-1];
s3.local = [1 0 0;0 1 0]; %coordinate frame of pixels
s3.pixelSize = 0.01;

% trace 1e3 times in 2 seconds
source = sourceColumn([0,0,0],[0,0,1],[1,0,0],10,40);
surfaces = {s1, s2, s3};

tic
trace = raytrace(source, surfaces);
toc
%
clf;
subplot(2,2,[1 2]);
plotSurfaces(trace);
axis equal; sideview
hold on;plotRays(trace,'b');hold off;

%%
subplot(2,2,3);
[wfe,mask] = pupilWFE(trace{end});
displayPupil(wfe,mask);colorbar
title(sprintf('RMS wfe %g \\mum',1e3*std(wfe(mask(:)))));
subplot(2,2,4);
oversample = 4;
psf = pupilPSF(wfe,mask,(500:10:600)*1e-6,[],100,20 * oversample,source.samplingDistance,.001/oversample);
psf = psf ./ max(psf(:));
imagesc(psf);
set(gca, 'XTick', oversample*(1:20), 'XTickLabel', 1:20)
set(gca, 'YTick', oversample*(1:20), 'YTickLabel', 1:20)
axis image; grid on

%% compare
% std(rays.position - finish.position,1)
% std(rays.direction - finish.direction,1)
