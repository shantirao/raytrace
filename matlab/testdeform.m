% test grid deformation

% primary mirror
PM=struct();
PM.label = 'PM';
PM.position = [0,0,0];
PM.direction= [0,0,1];
PM.cuy = 0; %flat
PM.K = 0;
PM.type = 1; % 0=reference, 1=reflect, 2=refract
PM.diameter = 100; %1m
PM.local = [1 0 0; 0 1 0];
PM.n = 1.0;

screen = struct();
screen.position=[0 0 100];
screen.direction=[0 0 -1];
screen.cuy = 0; %flat
screen.K = 0;
screen.type = 0; % 0=reference, 1=reflect, 2=refract
screen.local=[0 1 0; 1 0 0];
source = columnSource([0,0,100], [0,0,-1], [1,0,0], 50, 25);

%% before spot diagram
figure(1); clf;
[rays1]=raytrace(source,PM);
plotSpot(rays1);

%% generate deformation
step = 25;
[u,v] = meshgrid(-50:step:50,-50:step:50);

%tilt
Z = -abs(u)/10;
PM.deformation = deformSurfnodes = dlmread('d:/proj/itb/models/m1gradient/nodes.txt');
sag = dlmread('d:/proj/itb/models/m1gradient/sag.txt');

m1nodes = nodes(sag(:,1),:);
figure(4);
M1d = M1;

M1d.deformation =ace(PM,[u(:) v(:)],[Z(:) 0*v(:) Z(:)]);
% PM.deformation = deformSurface(PM,[u(:) v(:)],[0*u(:) 0*v(:) abs(u(:).*v(:))/100];

figure(2); clf; mesh(u,v,Z); 
hold on; 
triplot(PM.deformation.triangulation);  
for i=1:size(PM.deformation.triangulation,1);
    a = PM.deformation.points(PM.deformation.triangulation(i,:),:);
    a = mean(a,1);
    d = PM.deformation.curl(i,:);
    l = [a;a+100*d];
    if d(2) < 0
        c = 'b';
    else
        c = 'r';
    end
    plot3(l(:,1),l(:,2),l(:,3),c);
end
hold off;

%% Test the raytracer
figure(3); clf;
source3 = struct;
source3.chief=1;
source3.opl=[0;0];
source3.position=[1 1 100; -1 -1 100];
source3.direction=[0 0 -1; 0 0 -1];
source3.n2=1;
source3.N=2;
source3.map=[1 0; 0 1];
source3.aperture=1;

[rays3]=raytrace(source3,PM);
plotSpot(rays3);
%% Imaging system summary
figure(4); clf;

[rays4]=raytrace(source,{PM,screen});
plotSpot(rays4);
% subplot(2,2,[1 2])
% plotSideView(trace,[0 0 1;0 1 0],'b');
% hold on;
% plotOpticsSideView({trace{2:end}},[0 0 1;0 1 0],'r');
% hold off;
% 
% subplot(2,2,3)
% plotSpot(rays2);
% 
% subplot(2,2,4)
% [pupil,mask,rmswfe] = displayWFE(rays2);
% title(['Nominal WFE: ' num2str(1e6*rmswfe) 'nm rms']);
% psf = displayPSF(pupil, mask);
% colorbar;

% rays.N         [1x1] number of rays, N
% rays.n2        [1x1] (current index of refraction)^2
% rays.position  [Nx3] ray origins
% rays.direction [Nx3] ray directions
% rays.chief     [1x1] (set by source creator) index of chief ray
% rays.opl       [Nx1] (set by intersection creator) current path length