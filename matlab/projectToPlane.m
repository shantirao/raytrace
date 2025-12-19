function [points,distance] = projectToPlane(vertex,direction,xyz)
distance = (xyz - vertex) * direction';
points = xyz - distance * direction;
end
