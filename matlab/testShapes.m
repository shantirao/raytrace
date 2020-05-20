%% test pie
m1 = prescription{1};

pie.type = 'pie';
pie.radius = [1000 2500];
pie.edges = [2 0; 1 1.9];
pie.origin = m1.center;
pie.gap = [100 200];

source = columnSource(aperture,99,1);
options.negative=true;

figure(1); clf;

m1.aperture = [];
rays=raytrace(source,m1,options);
plotSpot(rays,10,'b'); 

m1.aperture = pie;
prays=raytrace(source,m1,options);

hold on; 
plotSpot(prays,15,'r','filled'); 

plotVector(surfaces{1}.center(:,1:2),[1e3 0],'b');
plotVector(surfaces{1}.center(:,1:2),[0 1e3],'g');
hold off;

%%
figure(2);clf;
plotApertures(surfaces{1});axis equal;

%% scratchpad
% label points
% hold on;arrayfun(@(i)text(uv(i,1),uv(i,2),num2str(i)),1:length(uv));hold off;

m1.segments = makePieSegments(m1,[50 2000 1000],50,6,2*pi/8);

% surface center needs to be offset by the parent offset or something
rays2=raytrace(source,m1,options);
plotSpot(rays2,20);
