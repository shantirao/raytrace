function [ pupil,mask] = pupilWFE(rays, rmttp)
%pupilWFE pupilWFE(rays) returns a pupil-referred optical path difference
% rmttp is whether to remove the Tip/Tilt/Piston (defaults to false)
% false means it removes the chief ray OPL only.
%map
% 
% opl = rays.opl(rays.valid);
% if isfield(rays,'chief') && rays.valid(rays.chief)
%     opd = opl - opl(rays.chief);
% else
%     opd = opl - mean(opl);
% end
% wfe = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),opd));
% mask = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),true));
% 
% % not the most elegant way of doing it
% if nargin > 1 && rmttp
%     pupil = rmTTP(wfe,mask);
% else
%     pupil = wfe;
% end

[pupil,mask] = pupilOPL(rays,rmttp);

if ~rmttp
  if isfield(rays,'chief') && rays.valid(rays.chief)
    chiefOPL = rays.opl(rays.chief);
  else
    chiefOPL = mean(pupil(mask(:)));
  end
  pupil = pupilOffset(pupil,mask,-chiefOPL);
end

end
