function [ru,zu,rv,zv] = surfaceCrossSection(surface)

aperture = max(abs(surface.aperture(:)));
N = 99;
M = floor(N/2);
R = (-M:M)'*aperture/M;

rays.chief = ceil(N/2);
rays.opl = zeros(N,1);
rays.position = bsxfun(@plus,surface.position+surface.direction*1000,bsxfun(@times,R,repmat(surface.local(1,:),N,1)));
rays.local = surface.local;
rays.direction = repmat(-surface.direction,N,1);
rays.n2 = 1;
rays.N = N;
rays.valid = true(N,1);

source1 = rays;
rays.position = bsxfun(@plus,surface.position+surface.direction*1000,bsxfun(@times,R,repmat(surface.local(2,:),N,1)));
source2 = rays;

rays1 = raytrace(source1,surface);
rays2 = raytrace(source2,surface);

%% elevation

Z1 = bsxfun(@minus,rays1.position,surface.position) * surface.direction';
Z2 = bsxfun(@minus,rays2.position,surface.position) * surface.direction';

ru = R(rays1.valid);
zu = Z1(rays1.valid);
rv = R(rays2.valid);
zv = Z2(rays2.valid);