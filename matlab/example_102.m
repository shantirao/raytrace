% Example 101: column source through a spherical biconvex lens

%% setup
% aperture
s0=struct();
s0.position = [0,0,0];
s0.direction = [0,0,1];
s0.aperture = 10;

% lens front surface
s1=struct();
s1.position = [0,0,40];
s1.direction=[0,0,1]; % points in the direction of curvature, and set convex=true
                      % (compatible with Zemax)
                      % or, points the other way, and set cuy = -1/50
                      % (compatible with MACOS)
s1.convex = true;
s1.cuy = 1/50;
s1.type = 'refract';
s1.index = 1.518522387620793;


% lens back surface -- it's concave, from the ray's point of view
s2=struct();
s2.position = [0, 0, 50];
s2.direction=-s1.direction;
s2.cuy = +1/50;
s2.type = 'refract';
s2.index = 1;

% Stop
s3 = findFocus(s0,{s1,s2});

%% trace 1e3 times in 2 seconds

source = sourceColumn(s0,5);
surfaces = {s1, s2, s3};

tic
trace = raytrace(source, surfaces);
toc

%%
clf;
plotSurfaces(trace);
axis equal; view(normr([1,0,1]));camroll(-90);
hold on;plotRays(trace,'b');hold off;
hold on; plotSurfaces(trace); hold off;
