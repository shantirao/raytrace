function scale = scaleUnits(units)
 %default to mm
if nargin < 1
    scale = 1e3;
elseif isnumeric(units)
    scale = 1/units;
elseif strcmp(units,'mm')
    scale =  1e3;
elseif strcmp(units,'um')
    scale =  1e6;
elseif strcmp(units,'nm')
    scale =  1e9;       
elseif strcmp(units,'pm')
    scale =  1e12;
elseif strcmp(units,'A')
    scale =  1e10;
else
    scale = 1;
end
end