function [ rays ] = sourceColumn(position, direction, x, radius, Nr, RefIndex)
%diverging Generates a parallel, cylindrical ray bundle with uniform spacing 
% Short syntax: columnSource(surface, Nrays, [RefIndex])
% surface        optical surface with position, direction, aperture, and local X axis
% Nr             [1x1] number of rays
% RefIndex       [1x1] scaling for ray direction (optional)
% Long syntax: columnSource(position, direction, x, radius, Nr, [RefIndex])
% position       [1x3] center of chief ray
% direction      [1x3] direction of chief ray
% x              [1x3] local x axis
% radius         [1x1] radius of pupil
% Nr             [1x1] integer number of rays on radius
% RefIndex       [1x1] scaling for ray direction
% returns rays: structure
% rays.N         [1x1] number of rays, N
% rays.n2        [1x1] (current index of refraction)^2
% rays.position  [Nx3] ray origins
% rays.direction [Nx3] ray directions
% rays.chief     [1x1] (set by source creator) index of chief ray
% rays.opl       [Nx1] (set by intersection creator) current path length

rays = struct;

if isstruct(position) % surface, Nrays, [RefIndex]
    srf = position;    
    RefIndex = 1;
    if nargin < 2
        Nr = 99;
    else
        Nr = direction; %num rays        
        if nargin > 2 
            RefIndex = x;
        end
    end
    if isfield(position,'units')
        rays.units = position.units;
        rays.display = position.units;
    end
    if isfield(position,'display')
        rays.display = position.display;
    end    
    if isfield(position,'index')
        RefIndex = position.index;
    end
    position = srf.position;
    direction = srf.direction;
    if isfield(srf,'local')
        x = srf.local(1,:);
    else
        x = [1 0 0];
    end
    if isfield(srf,'aperture')
        radius = max(srf.aperture);
    else
        radius = srf.diameter/2;
    end
else
    if nargin < 6
        RefIndex = 1;
    end
end

direction = normr(direction);
local = surfaceLocal(direction);
x = local(1,:);
y = local(2,:);
% y = normr(cross(direction,x));

if Nr == 0
    rays.chief = 1;
    rays.opl = 0;
    rays.position = position;
    rays.local = local;
    rays.direction = direction * RefIndex;
    rays.n2 = RefIndex^2;
    rays.N = 1;
    rays.map = [0 0];
    rays.aperture = radius;
    rays.valid = true;
else
    % uniform grid sampling by angle
    dU = radius*(-Nr:Nr)/Nr;
    rays.samplingDistance = dU(2)-dU(1);
    [Ux,Uy] = meshgrid(dU,dU);
    include = (Ux.^2+Uy.^2) <= radius^2;
    [r,c] =find( include );

    M = numel(Ux); % (2*Nr+1)^2
    N = nnz(include);
    Ux = Ux(include);
    Uy = Uy(include);

    rays.chief = ceil((2*Nr+2)*Nr+1);
%     rays.opl = zeros(N,1);
    rays.position = bsxfun(@plus,position,Ux*x + Uy*y);
    rays.local = local;
    rays.direction = repmat(direction * RefIndex,N,1);
    rays.n2 = RefIndex^2;
    rays.N = N;
    rays.map = [r c];
    rays.maskSize = [max(r), max(c)];
    % rays.mask = sparse(r,c,1);
    rays.aperture = radius;
    rays.valid = true(N,1);
end
end

