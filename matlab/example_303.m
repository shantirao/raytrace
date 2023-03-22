% Example 303: asphere, facing the other way.

% [start, finish] = importVV('example_301_VV.txt');

%% setup
% lens front surface
s1=struct();
s1.position = [0,0,10.427852];
s1.direction=[0,0,-1];
s1.curvature = -1/50; % negative means convex, or use the outside of the surface
s1.K = -2.30591024;
s1.type = 2; % 0=reference, 1=reflect, 2=refract
s1.asphere = [1e-8 -2e-13 3e-17 -4e-21];
s1.index = 1.518522387620793; % bk7

% lens back surface, flat
s2=struct();
s2.position = s1.position + [0,0,10];
s2.direction=[0,0,-1];
s2.type = 2;
s2.index = 1;

% Stop
s3=struct();
s3.position = s2.position+ [0,0,20];
s3.direction=[0,0,-1];
s3.type = 0;


%% trace 1e3 times in 2 seconds
source = sourceColumn([0,0,0],[0,0,1],[1,0,0],10,10);
surfaces = {s1, s2, s3};

tic
trace = raytrace(source, surfaces);
toc
%%
clf;
plotSurfaces(trace);
axis equal; view(normr([1,0,1]));camroll(-90);
hold on;plotRays(trace,'b');hold off;

%% compare
% std(rays.position - finish.position,1)
% std(rays.direction - finish.direction,1)
