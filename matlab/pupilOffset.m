function pupil = pupilOffset(pupil,mask,offset,fn,unmasked)
if numel(offset) == numel(pupil)
    if nargin > 3
        pupil(mask(:)) = fn(pupil(mask(:)),offset(mask(:)));
    else
        pupil(mask(:)) = pupil(mask(:)) + offset(mask(:));
    end
else
    if nargin > 3
        pupil(mask(:)) = fn(pupil(mask(:)),offset);
    else
        pupil(mask(:)) = pupil(mask(:)) + offset;
    end
end
if nargin> 4
    pupil(~mask(:))=unmasked;
end
end