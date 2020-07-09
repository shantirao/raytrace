function segments = makeHexSegments(parent,center,sizePP,gap,rotation,numRings,style)
% usage: m1.segments = makeHexSegments(m1,[0,0],500,50,0,2)
% works for up to 19 hex segments. After that, YMMV
% size is point-to-point size of the segments
s0 = struct(); %s0 becomes the prototype
for f = fieldnames(parent)' %have to use the transpose thingy.
    if ~strcmp(f{1},'segments')
        s0.(f{1}) = parent.(f{1});
    end
end
if isfield(parent,'center')
    center = parent.center(1:2) + center;
end
s0.center = center;

if nargin > 6
    s0.style = style;
end

name = '';
if isfield(parent,'name'), name = parent.name; end
wedge = 2*pi/6;


% h is height above 
[h, z] = elevation(s0);
z = -z; %segment Z is away from COC
s0.center(3) = h;
s0.aperture = makeHexagon(rotation,sizePP/2);
y = normr(surfaceLocalToGlobal(s0,s0.local(end,:)));
x = normr(cross(y,z));
y = normr(cross(z,x));
s0.reference = [x; y; z]; %= [parent.local;parent.direction];
s0.name = [name '.A0'];

spacing = sizePP * cos(wedge/2) + gap;
segments = {s0};


if numRings > 1
    centers = makeHexagon(rotation+wedge/2, spacing);
    for i=1:6
        s = s0;
        s.center = center+centers(i,:) ;
        [h, z] = elevation(s);
        z = -z;
        s.center(3) = h;
        s.aperture = makeHexagon(rotation,sizePP/2);
        s.name = [name '.A' num2str(i)];        
        y = normr(surfaceLocalToGlobal(s0,s0.center) - surfaceLocalToGlobal(s,s.center));
        % center is is in local coordinates, reference is is global
        x = normr(cross(y,z));
        y = normr(cross(z,x)); % make Y be in-plane     
%         z = normr(cross(x,y));
        s.reference = [x; y; z];
        segments{end+1} = s;
    end
end

if numRings > 2
    centers = makeHexagon(rotation,2* spacing * cos(wedge/2));
    for i=1:6
        s = s0;
        s.center = center+centers(i,:);
        [h, z] = elevation(s);
        z = -z;
        s.center(3) = h;
        s.aperture = makeHexagon(rotation,sizePP/2);
        s.name = [name '.B' num2str(i)];
        y = normr(surfaceLocalToGlobal(s0,s0.center) - surfaceLocalToGlobal(s,s.center));
        x = normr(cross(y,z));
        y = normr(cross(z,x));        
        s.reference = [x; y; z];
        segments{end+1} = s;
    end
    centers = makeHexagon(rotation+wedge/2,2* spacing );
    for i=1:6
        s = s0;
        s.center = center+centers(i,:);
        [h, z] = elevation(s);
        z = -z;
        s.center(3) = h;
        s.aperture = makeHexagon(rotation,sizePP/2);
        s.name = [name '.C' num2str(i)];
        y = normr(surfaceLocalToGlobal(s0,s0.center) - surfaceLocalToGlobal(s,s.center));
        x = normr(cross(y,z));
        y = normr(cross(z,x));        
        s.reference = [x; y; z];
        segments{end+1} = s;
    end
end

% for i=1:numel(segments)
%     segments{i}.center(3) = elevation(parent,segments{i}.center);
% end

function corners = makeHexagon(rotation,radius)
    theta = 2*pi/6*(0:5)+rotation;
    corners = radius*[cos(theta') sin(theta')];
end

end

