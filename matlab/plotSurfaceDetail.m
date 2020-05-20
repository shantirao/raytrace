function [ru,zu,rv,zv] = plotSurfaceDetail(surface)

[ru,zu,rv,zv] = surfaceCrossSection(surface);

% plot(R(rays1.valid),Z1(rays1.valid),R(rays2.valid),Z2(rays2.valid));
plot(ru,zu,rv,zv);
axis equal;
