function [focus] = findFocus(src,surfaces)
%% find the plane where parallel rays intersect
% inputs: surfaces is a cell array of optical surfaces, sourceAperture is a
% structure describing the ray origins
% outputs: pupil is the surface, index is where it goes in the surface 
% list, h is a plot handle. 

%% setup 
% src = sourceColumn(sourceAperture,0,1); % single ray at center of aperture
options.segments = false;
options.aperture = false;
%% find an image by tracing from the source aperture. 

ri = sourceColumn(src,9);
% ri = struct;
% ri.N = 2;
% ri.direction = repmat(src.direction,2,1);
% local = surfaceLocal(src);
% ri.position = [src.position;src.position+max(src.aperture/2)*local(1,:)]; % 1mm offset
% ri.valid = [true;true];
% ri.opl = [0;0];
% ri.chief =1;

r = raytrace(ri,surfaces,options);
r = r{end};
focus.name = 'Focus';
focus.position = lineIntersection(r.position,r.direction);
focus.direction = -r.direction(r.chief,:);
focus.local = surfaceLocal(focus);
end