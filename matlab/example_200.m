% Example 200: conic

% [start, finish] = importVV('example_200_VV.txt');

%% setup
% lens front surface
s1=struct();
s1.position = [0,0,10.427852];
s1.direction=[0,0,-1];
s1.cuy = -1/50;
s1.K = -2.30591024;
s1.type = 2; % 0=reference, 1=reflect, 2=refract
s1.n = 1.518522387620793;


% lens back surface, flat
s2=struct();
s2.position = s1.position + [0,0,10];
s2.direction=[0,0,-1];
s2.type = 2;
s2.n = 1;

% Stop
s3=struct();
s3.position = s2.position;
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