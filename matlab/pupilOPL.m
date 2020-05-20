function [ pupil,mask,tip,tilt] = pupilOPL(rays,rmttp)
%pupilOPL pupilOPL(rays) returns cumulative optical path difference
% rmttp is whether to remove the Tip/Tilt/Piston (defaults to false)

% pupil = full(sparse(rays.map(:,1),rays.map(:,2),rays.opl(:)));
% mask = full(sparse(rays.map(:,1),rays.map(:,2),true));
pupil = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),rays.opl(rays.valid),rays.maskSize(1),rays.maskSize(2)));
mask = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),true,rays.maskSize(1),rays.maskSize(2)));

sz = size(pupil);
if rays.maskSize(2) ~= sz(2)
    pupil = [pupil zeros(sz(1),rays.maskSize(2) - sz(2))];
    mask = [mask false(sz(1),rays.maskSize(2) - sz(2))];
end
sz = size(pupil);
if rays.maskSize(1) ~= sz(1)
    pupil = [pupil; zeros(rays.maskSize(1)-sz(1),sz(2))];
    mask = [mask; false(rays.maskSize(1)-sz(1),sz(2))];
end
% if isfield(rays,'surface') && isfield(rays.surface,'opdOffset')
%     pupil(mask) = pupil(mask) + rays.surface.opdOffset; 
% end
% not the most elegant way of doing it
if nargin > 1 && rmttp
    [pupil,tip,tilt] = rmTTP(pupil,mask);
else
    tip=[];tilt=[];
end

end
