function [pupil] = findVirtualPupil(surfaces,aperture,showFigures)
%% find an exit pupil that is an optical conjugate of the first mirror
% inputs: surfaces is a cell array of optical surfaces, aperture is a
% structure describing the ray origins (really only need one position and direction.)
% I wonder what happens when the virtual pupil wanders across the field
% outputs: pupil is the surface that you can trace backwards from the final surface

if nargin < 3
    showFigures = false;
end

%debug surfaces=telescope;
%% setup
src = sourceColumn(aperture,0); % single ray at center of aperture
options.segments = false;
options.aperture = false; % ignore holes
%% find an image by tracing from the source aperture.

% first, trace the source chief ray to the first surface
r = raytrace(src,surfaces{1},options); %t1{2}; % first surface intersection, instead of ;

if showFigures
    if isnumeric(showFigures)
        figure(showFigures); clf;
    elseif ~ishandle(gcf)
        figure();
    end

%    subplot(2,5,[1,6]);
    subplot(1,2,1);
    plotRays({src,r},'b');
    plotApertures(surfaces,true);
    xlabel("X"); ylabel("Y"); zlabel("Z");
    axis equal;
end

% add a bundle of rays that intersect the curved pupil surface, 100 urad
% offset in angle, and then trace those until the end

pos = r.position(r.chief,:);
dir = r.direction(r.chief,:);
local = surfaceLocal(surfaces{1});

rp = sourcePoint(pos, dir, local(1,:), 1e-5, 1, sqrt(r.n2));

%acos(dot(rp.direction(1,:) ,rp.direction(5,:)))
%acos(dot(rp.direction(2,:) ,rp.direction(4,:)))

% now trace that to the end
% re is where the chief ray strikes the focus plane
t2 = raytrace(rp,surfaces(2:end),options);
re = t2{end};

if showFigures
    hold on; plotRays(t2); hold off;
end

%should be 5 rays. Take two pairs and solve for the closest approach. It's probably a caustic, after all.
% 1,4  and 2,3
if showFigures
  subplot(1,2,2);
  for i=1:5
    plotLine(re.position(i,:),re.position(i,:)-re.direction(i,:)*10);
    if i==1,   hold on; end;
  endfor
  legend({"1","2","3","4","5"});xlabel("X"); ylabel("Y"); zlabel("Z");
  hold off;
  axis equal; grid on
end
%acosd(dot(re.direction(1,:) ,re.direction(5,:)))
%acosd(dot(re.direction(2,:) ,re.direction(4,:)))
%acosd(dot(re.direction(1,:) ,re.direction(2,:)))
%acosd(dot(re.direction(2,:) ,re.direction(4,:)))
%acosd(dot(re.direction(2,:) ,re.direction(3,:)))

[p1,d1] = lineLineClosest(re.position([1,3],:),re.direction([1,3],:));
[p2,d2] = lineLineClosest(re.position([2,3],:),re.direction([2,3],:));
[p4,d4] = lineLineClosest(re.position([4,3],:),re.direction([4,3],:));
[p5,d5] = lineLineClosest(re.position([5,3],:),re.direction([5,3],:));

d = min([d1(2),d2(2),d4(2),d5(2)]);

% pupil traces back along the chief ray -- have to generate a new pupil for each field point
pupil.name = 'pupil';
pupil.type = 'virtual';
pupil.radius = d;
pupil.position = re.position(re.chief,:) - d * re.direction(re.chief,:);
pupil.direction = re.direction(re.chief,:);

