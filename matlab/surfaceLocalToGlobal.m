function xyz = surfaceLocalToGlobal(s,points)
if nargin < 2
    points = [0,0];
end
    xyz = bsxfun(@plus,s.position,surfaceLocalToTangent(s,points));
end