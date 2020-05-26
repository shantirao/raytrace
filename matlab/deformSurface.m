function deformation = deformSurface(surface, points, displacement, varargin )
%deformSurface(surface, points, vectors) adds a freeform surface deformation
%Usage: surface.deformation = deform(surface, points, vectors)
%Input: surface [struct]: see raytrace.m
%       points [Nx2]: list of origin points, relative to the optic position
%                     call bsxfun(@minus,points,surface.position) to
%                     convert from global coordinates. If Nx3, will be
%                     projected onto the tangent plane
%       displacement [Nx3]: list of displacement vecotrs for those points
%       options:
%       'local': points are already in the local coordinate system; don't
%                translate. This is assumed if there are only 2 axes.
%
%Output: deformation.triangulation [Mx3]: delaunay triangulation, projected
%                                         onto the tangent plane
%        deformation.displacement [Nx1]: displacement, transformed into the
%                                        local coordinate system
%        deformation.points [Nx2]: original points (projected onto tangent
%        plane)
%        deformation.normal [Mx3]: facet normal vector
%        deformation.slope [Mx2]: facet slope 
%        deformation.vertexSlope [Nx2]: average slope at a vertex 

% removed:       deformation.curl [Nx3]: rotation vector

% have to do this in the local coodrinate system so mesh deformations
% survive 
translateCoordinates = (size(points,2) > 2);
for i=4:nargin
    if strcmp(varargin{i},'local')
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
    points = (points - surface.position) * surface.local';
end

DT = delaunayTriangulation(points);
d.points = points;
d.triangulation = DT;
d.displacement = displacement;

% compute a facet normal and curl vector for each triangle
N = size(DT,1); 

% d.curl = zeros(N,3);
d.normal = zeros(N,3);
d.slope = zeros(N,2);

% take advantage of Stokes' theorem: 
% integral((curl F) * dAarea) = integral(f * dLine)
% before is the triangle on the tangent plane
% after is the deformed triangle
% of course, the undistorted facet is always in the tangent plane so its
% normal is [0 0 1]. (duh!)
for i=1:N
    % easier -- before and after triangles define normals
    % cross product between those is the curl!
%     before = points(DT(i,:),:);
    after(:,1:2) = points(DT(i,:),:);
    after(:,3) = displacement(DT(i,:),:);
    
    %b = normr(cross(before(1,:)-before(2,:),before(1,:)-before(3,:)));
    a = normr(cross(after(1,:)-after(2,:),after(1,:)-after(3,:)));
    d.normal(i,:) = a;
%   curl = cross(b,a);
    d.slope(i,:) = - a(1:2)./a(3); % a is the normal. Slope of curve is -1/normal
end

M = size(points,1);
d.vertexSlope = zeros(M,2);

% given a point, which triangles is it in? average those slopes for all
% triangles intersecting at a vertex
for i=1:M
    touches = any(DT.ConnectivityList ==  i, 2); %which triangles include this point?
    % have ids of the triangles
    % now, want the angle-weighted mean of the slopes of those triangles
    a=[];
    b=[];
    triangles = DT.ConnectivityList(touches,:);
    for j = 1:size(triangles,1) % there's probably a more compact way to write this
        v = triangles(j,:); % 3 vertices, one of which is i
        if v(1) == i
            a(end+1,:) = points(v(2),:)-points(v(1),:);
            b(end+1,:) = points(v(3),:)-points(v(1),:);
        elseif v(2) == i
            a(end+1,:) = points(v(1),:)-points(v(2),:);
            b(end+1,:) = points(v(3),:)-points(v(2),:);
        else
            a(end+1,:) = points(v(1),:)-points(v(3),:);
            b(end+1,:) = points(v(2),:)-points(v(3),:);       
        end
    end
    % dot product to find angles
    c = acos(sum(a.*b,2)./sqrt(sum(a.*a,2) .* sum(b.*b,2)));
    s = c .* d.slope(touches,:);
    d.vertexSlope(i,:) = sum(s,1) / sum(c);
end

% and the edges also need a treatment
d.edges = edges(DT);
M = size(d.edges,1);
d.edgeSlope = zeros(M,2);
for i=1:M
    ID = edgeAttachments(DT,d.edges(i,1),d.edges(i,2));
    d.edgeSlope(i,:) = mean(d.slope(ID{:},:),1);
end

deformation = d;