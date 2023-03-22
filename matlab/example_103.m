% Example 103: column source through a spherical biconvex lens,
% plot polarization

%% setup
% aperture
s0=struct();
s0.position = [0,0,0];
s0.direction = [0,0,1];
s0.aperture = 10;
s0.local = [1,0,0]; % the initial polarization vector

% plane, at normal incidence
s1=struct();
s1.position = [0,0,40];
s1.direction= -s0.direction;
s1.type = 'refract';
s1.index = 1.518522387620793;

% concave
s2=struct();
s2.position = [0, 0, 50];
s2.direction = -s0.direction;
s2.cuy = +1/20;
s2.type = 'refract';
s2.index = 1;

% Stop
% s3 = findFocus(s0,{s1,s2})
s3=struct();
s3.position = [0,0,78];
s3.direction = -s0.direction;

%% trace 1e3 times in 2 seconds

source = sourceColumn(s0,5);
surfaces = {s1, s2, s3};

tic
trace = raytrace(source, surfaces);
clf; plotRays(trace,'b'); hold on; plotSurfaces(trace); hold off;
toc

%% polarization rotation
%clf
%rays = trace{end};
%scale = 1/2;
%xyz = rays.position;
%scatter(xyz(:,1),xyz(:,2),'b'); axis equal; grid on;
%hold on;
%for i=1:rays.N
%    xyz = rays.position(i,:);
%    xyz(2,:) = xyz + scale*rays.local(i,:);
%    plot(xyz(:,1),xyz(:,2),'r');
%end
%hold off;
%title('polarization rotation near focus')
