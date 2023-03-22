% Example 003: chromatic ray through source through a dispersive prism

% [start, finish] = importVV('example_001_VV.txt');

%% setup
% prism front surface
s1=struct();
s1.position = [0,20,20];
s1.direction=normr([0,1,-1]);
s1.type = 'refract';
s1.material = makeGlass('schott','F9');

% BS back surface
s2=struct();
s2.position = [0,20,20];
s2.direction=-normr([0,.3,1]);
s2.type = 'refract';
s2.index = 1;

% Stop
s3=struct();
s3.position = [0,0,40];
s3.direction=[0,0,-1];

%% trace 1e3 times in 2 seconds
N = 11;
source = struct();
source.position=repmat([0,10,0],N,1);
source.direction=repmat(normr([0,.4,1]),N,1);
source.valid=true(N,1);
source.wavelength=.000400+.0000300*(1:N)';
source.units='mm';

surfaces = {s1, s2, s3};

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
plotLine(s2.position*viewaxes,(s2.position + 10*cross(s2.direction,[1,0,0]))*viewaxes,'k');
axis image; xlabel('z'); ylabel('y');
hold off;
%% compare
% disp('accuracy')
% std(rays.position - finish.position,1)
% std(rays.direction - finish.direction,1)
