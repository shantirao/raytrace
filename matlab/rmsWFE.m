function [rmswfe, piston] = rmsWFE(pupil,mask) 
%     opd = pupil(:);
%     opd = opd(mask(:)); %reshape(pupil.*mask,1,numel(mask));

    % opd = opd(rays.valid);
    x = pupil(mask(:));
    rmswfe = std(x); %std(opd(find(mask)));
    if nargout > 1
        piston = mean(x);
    end
end

