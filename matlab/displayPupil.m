function [rmswfe,units] = displayPupil(pupil,mask,rays,clim)

if nargin > 2
    [scale, units] = displayScaleFactor(rays);
else
    [scale, units] = displayScaleFactor();
end
rmswfe = rmsWFE(scale*pupil,mask);
pupil(~mask(:))=NaN; % alternative to AlphaData for Octave

if nargin > 3
    img = imagesc(pupil,scale*clim);
else
    img = imagesc(pupil);
end
%set(img,'AlphaData',mask);
axis xy; axis equal;
set(gca,'XTick',[]);
set(gca,'YTick',[]);

if nargin > 2
    if isfield('label',rays)
        caption = rays.label;
    elseif isfield(rays,'surface') && isfield(rays.surface,'name')
        caption = rays.surface.name;
    else
        caption = '';
    end

    title(sprintf('%s %5.5g %s rms',caption,rmswfe,units));
end
