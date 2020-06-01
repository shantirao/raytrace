% Example 101: point source through a spherical biconvex lens

% [start, finish] = importVV('example_101_VV.txt');

%% setup
% lens front surface
s1=struct();
s1.position = [0,0,40];
s1.direction=[0,0,1]; % points in the direction of curvature, and set convex=true
                      % (compatible with Zemax)
                      % or, points the other way, and set cuy = -1/50
                      % (compatible with MACOS)
s1.cuy = 1/50;
s1.convex = true;
s1.type = 2; % 0=reference, 1=reflect, 2=refract
s1.n = 1.518522387620793;


% lens back surface -- it's concave, from the ray's point of view
s2=struct();
s2.position = [0, 0, 50];
s2.direction=-s1.direction;
s2.cuy = +1/50;
s2.type = 2;
s2.n = 1;

% Stop
s3=struct();
s3.position = [0,0,60];
s3.direction=-s2.direction; % flat doesn't care whether the direction points toward or away from the oncoming rays
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