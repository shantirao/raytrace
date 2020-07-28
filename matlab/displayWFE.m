function [pupil,mask,rmswfe] = displayWFE(rays, caption, rmttp, clim)
%DisplayWFE Helper function for plotting wavefront error maps
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
    [pupil1,mask1] = pupilWFE(rays{1}, rmttp);
    [pupil2,mask2] = pupilWFE(rays{2}, rmttp);
    pupil = pupil1 - pupil2;
    mask = mask1 & mask2;
    rays = rays{1}; % for printing units below
else
    [pupil,mask] = pupilWFE(rays, rmttp);
end

if nargin > 3
   [rmswfe,units] = displayPupil(pupil,mask,rays,clim);
else
   [rmswfe,units] = displayPupil(pupil,mask,rays);
end

if ~isempty(caption)
    title(sprintf('%s %5.5g %s rms',caption,rmswfe,units));
end
% title(sprintf('%s%5.5g %s rms',caption,scale*rmswfe,units));

% else
%     title([caption ': ' num2str(rmswfe) ' rms']);
% end
% colorbar; 
end

