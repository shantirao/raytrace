function [uv, e] = projectToTangentPlane(surface,xyz)

vertex = surface.position;
direction = surface.direction;
local = surface.local;

e = (xyz - vertex) * direction';
p = xyz - e * direction;

uv = p * local';

end