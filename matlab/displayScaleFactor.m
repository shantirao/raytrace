function [scale, units] = displayScaleFactor(rays)

if nargin == 0
    scale = 1;
    units = 'mm';
else
if iscell(rays)
    rays = rays{1}; 
end
scale = 1;
units = 'm';

% convert to meters
if isfield(rays,'units') 
    if strcmp(rays.units,'mm')
        scale = scale * 1e-3;
    elseif strcmp(rays.units,'um')
        scale = scale * 1e-6;
    elseif strcmp(rays.units,'nm')
        scale = scale * 1e-9;       
    elseif strcmp(rays.units,'pm')
        scale = scale * 1e-12;
    elseif strcmp(rays.units,'A')
        scale = scale * 1e-10;
    end
end

% convert to something other than meters
if isfield(rays,'display') 
    units = rays.display;
    if strcmp(units,'mm')
        scale = scale * 1e3;
    elseif strcmp(units,'um')
        scale = scale * 1e6;
    elseif strcmp(units,'nm')
        scale = scale * 1e9;       
    elseif strcmp(units,'pm')
        scale = scale * 1e12;
    elseif strcmp(rays.units,'A')
        scale = scale * 1e12;
    end
end
end