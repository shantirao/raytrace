function [pupil,mask,rmswfe] = plotRaytraceSummary2(trace,D)
subplot(2,3,[1 2 3])
rays = trace{end};
plotSideView(trace,[0 0 1;0 1 0],'b');
if isfield(rays,'NA')
title(sprintf('Numerical aperture %f, Diameter %f, EFL %f',rays.NA,D,D/rays.NA/2));
end
hold on;
plotOpticsSideView({trace{2:end}},[0 0 1;0 1 0],'r');
hold off;

subplot(2,3,4)
plotSpot(rays);
title(['Spot size: ' num2str(spotSize(rays)) 'mm rms']);

subplot(2,3,5)
[pupil,mask,rmswfe] = displayWFE(rays);
colorbar;

subplot(2,3,6)
psf = displayPSF(pupil, mask);
title(['Nominal WFE: ' num2str(1e6*rmswfe) 'nm rms']);
colorbar;