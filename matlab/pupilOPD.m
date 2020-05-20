function [ pupil,mask,rmswfe ] = pupilOPD(rays)
%pupilOPD pupilOPD(rays) returns a pupil-referred optical path difference
%map

opl = rays.opl(rays.valid);
opd = opl - mean(opl);
pupil = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),opd));
mask = full(sparse(rays.map(rays.valid,1),rays.map(rays.valid,2),1));
rmswfe = std(opd);

end

