function [ rays ] = sourceFan(position, direction, x, NA, nRays, RefIndex)
%diverging Generates a fat flan of 2n+1 rays with half-angle NA uniformly
%distributed in angle
% pointSource(position,direction,x,NA,Nrays) or
% pointSource(aperture,Nrays,RefIndex)
% aperture has position, direction, NA, local, n=RefIndex
% position       [1x3] center of central ray
% direction      [1x3] direction of central ray
% x              [1x3] local x axis of fan
% NA             [1x1] half-angle of cone
% nRays          [1x1] integer number of rays on radius
% RefIndex       [1x1] scaling for ray direction
% rays: structure
% rays.N         [1x1] number of rays, N
% rays.n2        [1x1] (current index of refraction)^2
% rays.position  [Nx3] ray origins
% rays.direction [Nx3] ray directions
% rays.central     [1x1] (set by source creator) index of central ray
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
c = -nRays:nRays;
Ux = U*(c)/nRays;

N = numel(Ux); 

rays.central = ceil((2*nRays+2)*nRays+1);
rays.position = repmat(position,N,1);
rays.direction = zeros(N,3);
rays.n2 = RefIndex^2;
rays.N = N;
rays.map = [ones(N,1),c'];
rays.maskSize = [N,1];
rays.local = [x;y];
rays.NA = NA;

rays.surface = struct;
rays.surface.position = position;
rays.surface.direction = direction;
rays.surface.local = rays.local;

z = direction' * RefIndex;
for i=1:N
   rays.direction(i,:) = rotationMatrix(y,Ux(i)) * z;
end

end

