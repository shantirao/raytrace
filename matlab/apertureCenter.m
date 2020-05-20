function [ center] = apertureCenter(surface)
%apertureCenter returns the center of the aperture in global coordinates,
%but relative to the vertex position

if isfield(surface,'center')
    if numel(surface.center) == 2
        center = surface.center*surface.local;
    else
        center = surface.center(1:2)*surface.local + surface.center(3)*surface.direction;
    end            
else
    center = [0,0,0];
end

end

