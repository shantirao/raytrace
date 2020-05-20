function [pupil,tip,tilt,piston] = rmTTP( pupil, mask )
%rmTTP Removes TTP from a pupil mask
%   Calculates overlap integral with tip and tilt masks, then subtracts
%   them. Seemed faster to rewrite than to find the previous
%   implementation.
%   divide tip and tilt by the optic radius to get the angle in radians

s = size(pupil);

if nargin < 2
    mask = (pupil ~= 0);
end

piston = mean(pupil(mask(:)));

pupil = mask .*(pupil - piston);

m = s(2)-1;
x = 2*((0:m) - (m/2))/s(2);
tip = repmat(x,s(1),1);
n = s(1)-1;
x = 2*((0:n) - (n/2))/s(1);
tilt = repmat(x',1,s(2));

an = sum(sum(mask.*(tip.*tip)));

a = sum(sum(pupil .* tip));
% pupil(mask(:)) =  (pupil(mask(:)) - (a/an .* tip(mask(:)) ));
pupil = mask .* (pupil - (a/an .* tip ));

bn = sum(sum(mask.*(tilt.*tilt)));

b = sum(sum(pupil .* tilt));
% pupil(mask(:)) = ( pupil(mask(:)) - (b/bn .* tilt(mask(:))));
pupil = mask .* ( pupil - (b/bn .* tilt));

%pupil(mask(:)) = pupil(mask(:)) - piston);

tip = a/an;
tilt = b/bn;
end

