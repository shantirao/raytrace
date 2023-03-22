% Example 002: column source through a tilted rectangular prism

% [start, finish] = importVV('example_001_VV.txt');

%% setup
% BS front surface
s1=struct();
s1.position = [0,0,70];
s1.direction=normr([0,-1,1]);
s1.type = 2; % 0=reference, 1=reflect, 2=refract
s1.index = 1.518522387620793;


% BS back surface
s2=struct();
s2.position = [0, 0, 80];
s2.direction=s1.direction;
s2.type = 2;
s2.index = 1;

% Stop
s3=struct();
s3.position = [0,0,120];
s3.direction=[0,0,1];
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
% disp('accuracy')
% std(rays.position - finish.position,1)
% std(rays.direction - finish.direction,1)
