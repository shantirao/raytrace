function xyz = surfaceLocalToTangent(s,points)
    u = points(:,1);
    v = points(:,2);
    if size(points,2) > 2
        w = points(:,3);
    else
        w = zeros(size(u)); %0 * u;
    end
    if isfield(s,'center')
        c = s.center; % back to the tangent plane projection
    else
        c = [0,0]; 
    end
    %
% 	xyz = ([s.local; s.direction] * [u+c(1); v+c(2); w])';
	xyz = ([u+c(1), v+c(2), w] *[s.local; s.direction]);
end