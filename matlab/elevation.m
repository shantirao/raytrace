function [ e, normal, xyz ] = elevation( surface, uv )
%elevation Calculates the elevation of a surface at positions on the
%tangent plane (relative to the surface center position, in local coordinates)
% points [Nx2] : defaults to surface.center [0,0]
% returns heights above tangent plane, normal vectors (in global
% coordinates) and positions (in global coordinates)
  
if ~isfield(surface,'cuy') && ~isfield(surface,'curvature') && ~isfield(surface,'radius')
    e = 0;
    normal = surface.direction;
    c = [0,0];
    if isfield(surface,'center')
        c = c + surface.center(1:2);
    end
    xyz = bsxfun(@plus,c,uv)*surfaceLocal(surface);
else
    if nargin < 2
        uv = [0,0];
    end

    if isfield(surface,'cuy')
    	z = abs(1/surface.cuy); %from ROC
    elseif isfield(surface,'curvature')
        z = 1/surface.curvature;
    elseif isfield(surface,'radius')
        z = surface.radius;
    end
    if isfield(surface,'convex') && surface.convex
        dir = -surface.direction;
    else
        dir = surface.direction;
    end
    N = size(uv,1);
    c = [0,0];
    if isfield(surface,'center')
        c = c + surface.center(1:2);
    end

    x = bsxfun(@plus,c,uv)*surfaceLocal(surface);
    rays.position = bsxfun(@plus,surface.position+z*dir,x);
    rays.direction = repmat(-dir,N,1);

    % flip the surface around so that we can trace from the tangent plane.
    % remember to flip the Z of the normal vector for this to make sense.

        surface.aperture = [];
        surface.segments = {};
        [rays,~,~,e] = raytrace(rays,surface);

        if nargout > 1
            %reverse the normal vector 
            normal = rays.normal; % - 2*bsxfun(@times,surface.direction,rays.normal*surface.direction');
        end
        if nargout > 2
            xyz = rays.position;
        end
    end
    
end

