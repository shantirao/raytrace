function [pupil,mask,rmswfe,tip,tilt] = displayOPL(rays, caption, rmttp, clim)
%displayOPL Helper function for plotting wavefront error maps
% rays can be a cell array, in which case the difference in wfe is plotted

if nargin < 2
    if isfield(rays,'surface') && isfield(rays.surface,'name')
        caption = [rays.surface.name ' '];
    else
        caption = '';
    end
else
    caption = [caption ' '];
end
if nargin < 3
    rmttp = false;
end

if iscell(rays)
    [pupil1,mask1] = pupilOPL(rays{1}, rmttp);
    [pupil2,mask2] = pupilOPL(rays{2}, rmttp);
    pupil = pupil1 - pupil2;
    mask = mask1 & mask2;
    rays = rays{1}; % for printing units below
else
    [pupil,mask,tip,tilt] = pupilOPL(rays, rmttp);
end


[scale, units] = displayScaleFactor(rays);

rmswfe = rmsWFE(pupil,mask);

wfe = pupilOffset(pupil,mask,-mean(pupil(mask(:))));

if nargin > 3
    img = imagesc(scale*wfe,scale*clim); 
else
    img = imagesc(scale*wfe); 
end
set(img,'AlphaData',mask)
axis image; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);
% if isfield(rays,'units')

if rmttp
    title(sprintf('%s %.5g %s rms, tip %g, tilt %g ',caption,scale*rmswfe,units,tip/rays.aperture,tilt/rays.aperture));
else    
    title(sprintf('%s %.5g %s rms',caption,scale*rmswfe,units));
end
% else
%     title([caption ': ' num2str(rmswfe) ' rms']);
% end
% colorbar; 
end

