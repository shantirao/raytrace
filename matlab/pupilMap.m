function [ pupil,mask] = pupilMap(rays,setnan)

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

if nargin > 1 && setnan
  pupil(~mask(:)) = NaN;
end

end
