function [ e, normal, xyz ] = elevation( surface, uv )
%elevation Calculates the elevation of a surface at positions on the
%tangent plane (relative to the surface center position, in local coordinates)
% points [Nx2] : defaults to surface.center [0,0]
% returns heights above tangent plane, normal vectors (in global
% coordinates) and positions (in global coordinates)
if ~isfield(surface,'cuy')
    e = 0;
    normal = surface.direction;
    c = [0,0];
    if isfield(surface,'center')
        c = c + surface.center(1:2);
    end
    xyz = bsxfun(@plus,c,uv) * surface.local;
else
    if nargin < 2
        uv = [0,0];
    end

    z = abs(1/surface.cuy); %from ROC
    N = size(uv,1);
    c = [0,0];
    if isfield(surface,'center')
        c = c + surface.center(1:2);
    end

    x = bsxfun(@plus,c,uv) * surface.local;

    rays.position = bsxfun(@plus,surface.position+z*surface.direction,x);
    rays.direction = repmat(-surface.direction,N,1);

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

