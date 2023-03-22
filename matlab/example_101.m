% Example 101: point source through a spherical biconvex lens
% compare to a Zemax model calculated by Norbert Sigrist
% difference should be < 1e-14 mm in position, < 1e-15 in direction
% (only 15 digits are provided in the data file)

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
s2.cuy = 1/50;
s2.type = 'refract';
s2.index = 1;

% Stop
s4=struct();
s4.position = [0,0,60];
s4.direction=[0,0,1];
s4.type = 0;

%%
[start, finish] = importVV('example_101_VV.txt');
trace = raytrace(start,{s1,s2,s4});
rays = trace{end};

plotSurfaces(trace);
axis equal; view(normr([1,0,1]));camroll(-90);
hold on;plotRays(trace,'b');hold off;
hold on; plotSurfaces(trace(2:end)); hold off;

disp('mean error')
mean(abs(rays.position - finish.position),1)
mean(abs(rays.direction - finish.direction),1)
disp('range')
max(abs(rays.position - finish.position),[],1)
max(abs(rays.direction - finish.direction),[],1)
disp('std')
std(rays.position - finish.position,1)
std(rays.direction - finish.direction,1)

%%
clf;
plotSpot(trace);
