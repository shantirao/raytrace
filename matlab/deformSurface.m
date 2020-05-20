function deformation = deformSurface(surface, points, displacement, varargin )
%deformSurface(surface, points, vectors) adds a freeform surface deformation
%Usage: surface.deformation = deform(surface, points, vectors)
%Input: surface [struct]: see raytrace.m
%       points [Nx3]: list of origin points, relative to the optic position
%                     call bsxfun(@minus,points,surface.position) to
%                     convert from global coordinates. Doesn't need to be
%                     projected onto the tangent plane.
%       displacement [Nx3]: list of displacement vecotrs for those points
%       options:
%       'local': points are already in the local coordinate system; don't
%                translate. This is assumed if there are only 2 axes.
%
%Output: deformation.triangulation [Mx3]: delaunay triangulation, projected
%                                         onto the tangent plane
%        deformation.displacement [Nx3]: displacement, transformed into the
%                                        local coordinate system
%        deformation.tilt [Nx3]: rotation vector

% have to do this in the local coodrinate system so mesh deformations
% survive 
translateCoordinates = (size(points,2) > 2);
for i=4:nargin
    if ~strcmp(varargin{i},'local')
        translateCoordinates = false;
    end
end

if length(points) < 3
    error('deformSurface: %s','need more than 3 points')
end
    
if ~isfield(surface,'local') || numel(surface.local) ~= 6
    error('deformSurface: %s','surface must have a 2D local coordinate system')
end

d = struct();

if translateCoordinates
    points = bsxfun(@minus,points,surface.position);
    projection = points * surface.local';
    uvw = [surface.local; cross(surface.local(1,:),surface.local(2,:))];
    points = points * xyz';
else
    projection = points;
    if size(points,2) == 2
        points(:,3)= 0;
    end
end

DT = delaunayTriangulation(projection);
d.points = points;
d.triangulation = DT;
d.displacement = displacement;

% compute a curl vector for each triangle
N = size(DT,1);
d.curl = zeros(N,3);

% take advantage of Stokes' theorem: 
% integral((curl F) * dAarea) = integral(f * dLine)
%area = 
for i=1:N
    % easier -- before and after triangles define normals
    %cross product between those is the curl!
    before = points(DT(i,:),:);
    after = before + displacement(DT(i,:),:);
    
    b = cross(before(1,:)-before(2,:),before(1,:)-before(3,:));
    a = cross(after(1,:)-after(2,:),after(1,:)-after(3,:));
    
    d.curl(i,:) = cross(normr(b),normr(a));
    
%     forces = displacement(DT(i,:),:);
%     positions = forces + points(DT(i,:),:);
%    
%     a = positions(2,:)-positions(1,:);
%     b = positions(2,:)-positions(3,:);
%     %actually twice the line integral because the area is half the
%     %cross-product
%     curl = forces(1,:)+forces(2,:) / norm(a) + ...
%         forces(2,:)+forces(3,:) / norm(b) + ...
%         forces(3,:)+forces(1,:) / norm(positions(1,:)-positions(3,:));
%     
%     if length(a) == 3
%         area = norm(cross(a,b));
%     else
%         area = abs(a(1)*b(2) - a(2)*b(1));
%     end
%     
%     d.curl(i,:) = curl/area;
end

deformation = d;