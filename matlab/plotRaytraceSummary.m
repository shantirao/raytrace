function [pupil,mask,rmswfe] = plotRaytraceSummary(lightpath,D)
  if nargin < 2
    D = lightpath{1}.aperture*2;
  endif
  subplot(2,3,[1 2 3])
  rays = lightpath{end};
  plotSideView(lightpath,[0 0 1;0 1 0],'b');
  if isfield(rays,'NA')
  title(sprintf('Numerical aperture %f, Diameter %f, EFL %f',rays.NA,D,D/rays.NA/2));
  end
  hold on;
  plotOpticsSideView({lightpath{2:end}},[0 0 1;0 1 0],'r'); axis equal;
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
