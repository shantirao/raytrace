function [rmswfe, piston] = rmsWFE(pupil,mask) 
%     opd = pupil(:);
%     opd = opd(mask(:)); %reshape(pupil.*mask,1,numel(mask));

    % opd = opd(rays.valid);
if isstruct(pupil) % rays
    opd = pupil.opl(pupil.valid(:));
    piston = mean(opd);
    rmswfe = std(opd-piston,1);
else
    x = pupil(mask(:));
    rmswfe = std(x,1); %std(opd(find(mask)));
    if nargout > 1
        piston = mean(x);
    end
end

