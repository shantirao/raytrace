function [ points ] = transform3D( points,T )
%transform3D Apply a 4x4 rigid body transformation matrix to a set of points.

points = T * [points'; ones(1,size(points,1))];
points = points(1:3,:)';

end

