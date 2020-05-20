function [pupil,mask] = trimPupil(pupil,mask)

if nargin<2, mask = pupil; end;
i=1;
ii=size(pupil,1);
j=1;
jj=size(pupil,2);

while i < ii && ~any(mask(i,:)), i = i+1; end
while i < ii && ~any(mask(ii,:)), ii = ii-1; end

while j < jj && ~any(mask(:,j)), j = j+1; end
while j < jj && ~any(mask(:,jj)), jj = jj-1; end

pupil=pupil(i:ii,j:jj);
mask=boolean(mask(i:ii,j:jj));

return