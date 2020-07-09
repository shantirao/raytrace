function [trace, distance, projection, elevation] = raytrace(rays, surface, options)
%raytrace(rays,surfaces,options) Finds where rays hit a surface. 
% {rays} = raytrace2(start,{surface list},options)
%Input: rays [struct] 
%       rays.position [Nx3]: positions
%       rays.direction [Nx3]: direction * index
%       rays.n: index
%       rays.n2: index^2
%       rays.local [Nx3]: local horizontal polarization vector, or TE
%       	vector. Each interface involves a coordinate frame rotation
%       	depending on the plane of incidence
%       rays.polarization [Nx4]: Stokes parameters S0..3 for the rays 
%           Default [ones(N,1) zeros(N,3)] for unpolarized light
%
%       surface [struct] or cell array {struct, struct...}
%           if a cell array, the returned {rays} is a raytrace series
%       surface.type [1x1]: surface type
%                       0: stop or reference (don't change direction)
%                       1: reflective
%                       2: refractive
%                       4: polarization - reflection if surface.NK is
%                       defined; otherwise reflection and transmission.
%       surface.display [text]: override display units for this surface 
%                       nm, um, mm, m
%       surface.segments {}: if present, nonsequential raytrace of segments
%                        could mix reflective, refractive, and stops.
%                        Segments don't necessarily have the same
%                        prescription but the container's surface.center
%                        field propagates through to make it easy to move
%                        off-axis apertures around
%       surface.n [1x1]: next index of refraction
%       surface.NK [1xM]: complex mirror surface index for M wavelengths
%       surface.position [1x3]: vertex point
%       surface.direction [1x3]: surface normal at vertex
%       surface.cuy [1x1]: curvature in y (semi-latus rectum) = 1/ROC
%                          cuy > 0: direction points towards   curved surface
%                          cuy < 0: direction points away from curved surface
%                          cuy = 0: Plane
%       surface.curvature [1x1]: same as cuy
%       surface.radius [1x1]: alternative to curvature, radius = 1/cuy
%       surface.K [1x1]: conic constant
%                          K = 0: for a sphere
%       surface.center [1x2]: center of aperture on the tangent plane
%                             (defaults to 0,0)
%       surface.local [2x3]: local coordinate system [ x ; y ], both
%                            orthogonal to surface.direction
%       surface.aperture [1x1]: aperture radius 
%                        [1x2]: [inner, outer] for an annulus
%                        [1x3]: [inner, semiMinor, semiMajor] ellipse (not yet)
%                        [2x2]: [u1 v1; u2 v2] rectangle (local coords)
%                        [Mx2]: [u1 v1; u2 v2; ... ] polygon (local coords)
%       surface.aperture.type == 'pie':
%                               .radius = [inner outer]
%                               .edges = [x1 y1; x2 y2] unit vectors describing the edges 
%                               .gap = [g1, g2] edge offset
%                               .origin = [u1 v1] 
%                                   parent center, not the segment's center
%                                   of rotation
%       surface.asphere [1xN]: even asphere coefficients a2 s^2 + ...
%                                  a4 s^4 + a6 s^6 + ...
%       surface.zernike [1xN]: Zernike coefficients (not implemented)
%       surface.deformation [struct]:  see deformSurface.m 
%            deformation.origins [3xN]: mesh points that lie on surface
%            deformation.displacement[3xN]: vector displacement of each coordinate
%            deformation.displacement[1xN]: optical sag of each coordinate
%            deformation.projection  [2xN]: origins, projected onto the tangent plane
%            deformation.triangles   [3xM]: indices of triangles
%       surface.pixelSize [1x1]: if you want the rays to map to pixel
%                          coordinates
%       surface.transform [4x4]: coordinate transform in case the optic
%                                was rotated or moved and the mesh wasn't.
%           [x1;y1;z1;~]=[RotationMatrix [dx;dy;dz]; 0 0 0 1]*[x0;y0;z0;1]
%       surface.reference [3x3]; [x; y; z] axes of reference geometry in
%                                the global coordinate system, used for
%                                perturbations
%       surface.reflectance [struct]
%               reflectance.angle [2xM]: indices of triangles
%
%       options.aperture [bool]: [true] use aperture stops
%       options.negative [bool]: [false] include negative rays
%       options.segments [bool]: [true] include segments
%       options.debug [bool]: [true] include distance, projection,
%           elevation in return structure
%Output: If rays and surface are structures, then rays is a structure.
%        Otherwise, rays is a cell array of structures (or arrays of
%        structures).
%        rays [struct]
%        rays.position [Nx3]: RSI position
%        rays.distance [Nx3]: RSI
%        rays.opl [Nx1]: optical path length
%        rays.valid [Nx1]: boolean: valid ray or not?
%        rays.n2: index^2 after last surface
%        rays.n: index after last surface
%        rayys.pixels [Nx3]: which pixels the rays land on
%        rays.cos: rays.direction . surface.N / rays.n;
%        rays.status: >0: stopped at surface N
%                      0: invalid
%                     -1: no Ray-Surface-Intersection
%                     -2: Negative ray
%                     -3: Blocked (Aperture)
%
%        these are only evaluated for single optic traces
%        distance   [Nx3]: distance traveled/N
%        projection [Nx3]: projection of intersection point on tangent plane
%        elevation [Nx1]: distance to the tangent plane

%   to do:
%       rays.lambda [1xM]: optional vector of wavelengths to 
%       rays.phase [NxM]: optional complex amplitude & phase (wavelength-dependent)
%               initialize as 1+0i, multiply by -1 at each reflection
%       rays {cell array}: iterate trace over multiple ray bundles
% Add Zernike, mesh, grid perturbations
% GPU is very slow. parallel across CPUs is better

if nargin < 2
    error('Usage: raytrace({rays}, {surfaces})');
end

if nargin < 3
    options = struct;
    options.aperture = true;   % stop rays that fall outside the aperture mask
    options.negative = false;  % no negative rays
    options.segments = true;   % segmented optics
    options.parallel = false;  % parallel processing isn't always faster
    options.precision = 1e-14; % 1e-14 mm should be enough
    options.debug = false;
end

if iscell(rays) % parallel
    if options.parallel
        trace = parmap(@(r)raytrace(r,surface,options),rays);
    else
        trace = cell(numel(rays));
        for i=1:numel(surface)
        	trace{i} = raytrace(rays{i},surface,options);
        end
    end
elseif iscell(surface) % series
    N = numel(surface);
    trace = cell(1,N+1); 
    trace{1}=rays;
    source = rays;
    for i=1:numel(surface)
       rays = raytrace(rays,surface{i},options);
       if isfield(source,'display') && ~isfield(surface{i},'display')
            rays.display = source.display; %source units override surface units
       end
       trace{i+1} = rays;
    end
elseif (~isfield(options,'segments') || (isfield(options,'segments') && options.segments)) && ...
    isfield(surface,'segments') && iscell(surface.segments) && numel(surface.segments) 
% parallel
    if (isfield(options,'parallel') && options.parallel)
        [R, D, P, E] = parmap(@(s)raytrace(rays,s,options),surface.segments);
    else
        [R, D, P, E] = map(@(s)raytrace(rays,s,options),surface.segments);
    end
   
    if iscell(R)
        rays = R{1};
    else 
        rays = R;
    end
    rays.surface = surface;
    % use rays.valid to collapse the other traces
    
    distance = D{1};
    projection = P{1};
    elevation = E{1};
    rays.segment = 1*rays.valid;
    for i=2:numel(surface.segments) %combine in series
        r = R{i};
        d = D{i};
        p = P{i};
        e = E{i};
        index = find(r.valid); %slightly faster than a binary map
        distance(index,:) = d(index,:);
        projection(index,:) = p(index,:);
        elevation(index,:) = e(index,:);
        rays.opl(index,:) = r.opl(index,:);
        rays.position(index,:) = r.position(index,:);
        rays.direction(index,:) = r.direction(index,:);
        if isfield(r,'map'), rays.map(index,:) = r.map(index,:); end
        if isfield(r,'status'), rays.status(index,:) = r.status(index,:); end
        rays.segment(index,:) = i; 
        rays.valid(index,:) = true;
    end
    trace = rays;
else % propagate one ray bundle to one surface
    [trace, distance, projection, elevation] = propagateRays(rays, surface, options);
end

