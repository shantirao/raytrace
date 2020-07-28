function [pupil,newsurfaces,pupilIndex] = findExitPupil(sourceAperture,surfaces,display)
%% find an exit pupil that is an optical conjugate of the first mirror
% inputs: surfaces is a cell array of optical surfaces, sourceAperture is a
% structure describing the ray origins
% outputs: pupil is the surface, index is where it goes in the surface 
% list, h is a plot handle. 

if nargin < 3
    display = false;
end

%% setup 
src = sourceColumn(sourceAperture,0); % single ray at center of aperture
options.segments = false;
options.aperture = false;
%% find an image by tracing from the source aperture. 

% scale = src.aperture/100;

%% find a pupil
% look for the optical conjugate of the first surface before the last image
% plane

% first, trace the source chief ray to the first surface

r = raytrace(src,surfaces{1},options); %t1{2}; % first surface intersection, instead of ;

if display
    if isnumeric(display)
        figure(display); clf;
    elseif ~ishandle(gcf)
        figure();
    end
    
    subplot(2,5,[1,6]); 
    plotRays({src,r},'b');
    plotApertures(surfaces,true);
    axis equal;
end

% add a bundle of rays that intersect the curved pupil surface, 100 urad
% offset in angle, and then trace those until the end 

pos = r.position(r.chief,:);
dir = r.direction(r.chief,:);
local = surfaceLocal(surfaces{1});

rp = sourcePoint(pos, dir, local(1,:), 1e-4, 1, sqrt(r.n2));

% now trace that to the end 
t2 = raytrace(rp,surfaces(2:end),options);
r2 = t2{end};

if display
    hold on;
    plotRays(t2,'r');
end
%%
%pimg is where the chief ray strikes the focus plane

% now trace those rays backwards to find where they intersect
r3 = r2;
r3.opl = r3.opl*0;
r3.direction = -r2.direction;

for pupilIndex = numel(surfaces)-1:-1:2
    pLast = r3.position(r3.chief,:); 
    [ppupil,d] = lineIntersection(r3.position,r3.direction); %(1:2,:)
    if norm(ppupil - pLast) < norm(r3.position(r3.chief,:) - surfaces{pupilIndex}.position) 
        break % is the intersection between the current surface and the next surface?
    end
    r3 = raytrace(r3,surfaces{pupilIndex},options); % go on to the next surface
end

pupil = struct;
pupil.name = 'Pupil relay';
pupil.position = ppupil; % vertex tangent

%%

% [ppupil,d] = lineIntersection(r2.position,r2.direction);
disp(sprintf('Pupil at %d %d %d at position ',ppupil, pupilIndex+1));

if display
    scatter3(ppupil(1),ppupil(2),ppupil(3),'r');  
    text(ppupil(1),ppupil(2),ppupil(3),'Pupil relay','color','r','HorizontalAlignment','center'); 
    plot3([ppupil(1), ppupil(1); r2.position(1,1), r2.position(2,1)], ...
        [ppupil(2), ppupil(2); r2.position(1,2), r2.position(2,2)], ...
        [ppupil(3), ppupil(3); r2.position(1,3), r2.position(2,3)],'g');
    hold off;
end
%% new pupil surface

ROC = r3.opl(r3.chief) + norm(ppupil - r3.position(r3.chief,:)); %ppupil - pimg;

pupil.opdOffset = ROC; %dot(dV,direction); %negative means virtual
pupil.radius = -ROC; %sqrt(sum((dV).^2)); % positive = concave
pupil.direction = r3.direction(r3.chief,:); 
% pupil.convex = true;
pupil.local = surfaceLocal(pupil);
pupil.display = 'nm';

newsurfaces = { surfaces{1:pupilIndex}, pupil, surfaces{pupilIndex+1:end}};

pupilIndex = pupilIndex+1;
%%  plot summary
if display
    source = sourceColumn(sourceAperture,3,1);
    source.units = 'mm';
    source.display = 'nm';
    options.negative = true;
    trace = raytrace(source,surfaces,options);
    
    subplot(2,5,3); [p1,mask] = displayOPL(trace{end-1});
    
    subplot(2,5,4);	p2 = displayOPL(trace{end});
    
    subplot(2,5,5); 
    dp = p1-p2;
    dp(mask) = dp(mask) - dp(trace{end}.chief);
    img =  imagesc(1e6*dp); title(sprintf('Difference: %5.5g nm',1e6*std(dp(mask(:)))));
    set(img,'AlphaData',mask); axis image
 
    subplot(2,5,8); plotSpot(trace{end-2});  
    if isfield(surface,'name'),title(trace{end-2}.surface.name); end
    
    subplot(2,5,9); plotSpot(trace{end-1});  
    if isfield(surface,'name'), title(trace{end-1}.surface.name); end
    
    subplot(2,5,10); plotSpot(trace{end});  
    if isfield(surface,'name'),title(trace{end}.surface.name);end
    
    subplot(2,5,[2 7]); plotRays(trace,'b');plotSurfaces(trace);axis image;
    
    
end
end