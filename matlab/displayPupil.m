function rmswfe = displayPupil(pupil,mask,rays,clim)


rmswfe = rmsWFE(pupil,mask);

[scale, units] = displayScaleFactor(rays);
if nargin > 3
    img = imagesc(scale*pupil,scale*clim); 
else
    img = imagesc(scale*pupil); 
end
set(img,'AlphaData',mask)
axis image; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);
if isfield('label',rays)
    caption = rays.label;
else 
    caption = '';
end

title(sprintf('%s%5.5g %s rms',caption,scale*rmswfe,units));
