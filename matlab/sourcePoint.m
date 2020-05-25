function [ rays ] = sourcePoint(position, direction, x, NA, nRays, RefIndex)
%diverging Generates a cone of diverging rays with half-angle NA uniformly
%distributed in sine-angle
% sourcePoint(position,direction,x,NA,Nrays) or
% sourcePoint(aperture,Nrays,RefIndex)
% aperture has position, direction, NA, local, n=RefIndex
% position       [1x3] center of chief ray
% direction      [1x3] direction of chief ray
% x              [1x3] local x axis
% NA             [1x1] half-angle of cone
% nRays          [1x1] integer number of rays on radius
% RefIndex       [1x1] scaling for ray direction
% rays: structure
% rays.N         [1x1] number of rays, N
% rays.n2        [1x1] (current index of refraction)^2
% rays.position  [Nx3] ray origins
% rays.direction [Nx3] ray directions
% rays.chief     [1x1] (set by source creator) index of chief ray
% rays.opl       [Nx1] (set by intersection creator) current path length
if nargin < 6
    RefIndex = 1;
end
if isstruct(position)
    s = position;  
    if nargin>2, RefIndex=x; elseif isfield(s,'n'), RefIndex=s.n; end
    if nargin>1, nRays=direction; else nRays = 99; end
    position = s.position;
    direction = s.direction;
    x = s.local(1,:);
    if isfield(s,'NA'), NA = s.NA; else NA=pi/4; end
    if isfield(s,'units')
        rays.units = s.units;
        rays.display = s.units;
    end
    if isfield(s,'display')
        rays.display = s.display;
    end    
end

direction = normr(direction);
y = normr(cross(direction,x));

% uniform grid sampling by angle
U = NA/RefIndex;
dU = U*(-nRays:nRays)/nRays;

if U > pi/4
    error('pointSource only works up to 45º')
end

[Ux,Uy] = meshgrid(dU,dU);
include = (Ux.^2+Uy.^2) <= U^2;
[r,c] =find( include );

M = numel(Ux); % (2*Nr+1)^2
N = nnz(include);
Ux = asin(Ux(include));
Uy = asin(Uy(include));

rays.chief = ceil((2*nRays+2)*nRays+1);
rays.position = repmat(position,N,1);
rays.direction = zeros(N,3);
rays.n2 = RefIndex^2;
rays.N = N;
rays.map = [r c]; %full(sparse(r,c,1:N));
rays.maskSize = [max(r), max(c)];
rays.valid = true(N,1);

z = direction' * RefIndex;
for i=1:N
   th = -x*Uy(i) + y*Ux(i); % rotation vector
   angle = norm(th);
   axis = th/angle;
   rays.direction(i,:) = rotationMatrix(axis,angle) * z;
 %  rays.direction(i,:) = QRotMatrix(th) * z;
end
    
end

