function [ vectors ] = transformVector( vectors,T )
%transform3D Apply a 4x4 rigid body transformation matrix to directions.

vectors = T * [vectors'; zeros(1,size(vectors,1))];
vectors = vectors(1:3,:)';

end

