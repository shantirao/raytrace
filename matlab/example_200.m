% Example 200: conic

% [start, finish] = importVV('example_200_VV.txt');

%% setup - perfectly aligned conic lens has a very boring PSF
% lens front surface
s1=struct();
s1.position = [0,0,10.427852];
s1.direction=[0,0,-1];
s1.curvature = -1/50;
s1.convex = true;
s1.K = -2.30591024;
s1.type = 'refract';
s1.index = 1.518522387620793;


% lens back surface, flat
s2=struct();
s2.position = s1.position + [0,0,10];
s2.direction=[0,0,-1];
s2.type = 'refract';
s2.index = 1;

% Stop
s3=struct();
s3.position = s1.position + [0,0,100];
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
clf; subplot(2,2,[1 2]);
plotSurfaces(trace);
axis equal; sideview
hold on;plotRays(trace,'b');hold off;

%%
subplot(2,2,3);
[wfe,mask] = pupilWFE(trace{end});
displayPupil(wfe,mask);colorbar

subplot(2,2,4);
oversample = 5;
psf = pupilPSF(wfe,mask,(500:10:600)*1e-6,[],100,20*oversample,source.samplingDistance,s3.pixelSize/oversample);
psf = psf ./ max(psf(:));
imagesc(psf);
set(gca, 'XTick', oversample*(1:20), 'XTickLabel', 1:20)
set(gca, 'YTick', oversample*(1:20), 'YTickLabel', 1:20)
axis image; grid on
%% compare
% std(rays.position - finish.position,1)
% std(rays.direction - finish.direction,1)
