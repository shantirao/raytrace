function xyz = surfaceGlobalToTangent(surface, xyz)
% project a global point onto the surface's tangent plane
    p = surface.position - xyz;
    d = normr(surface.direction);
    l = (d * p')';
    delta = bsxfun(@times,d,l);
    xyz = bsxfun(@plus,xyz,delta);
end


