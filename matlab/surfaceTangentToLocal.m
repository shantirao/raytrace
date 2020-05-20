function uv = surfaceTangentToLocal(s,xyz)
% Like surfaceGlobalToLocal only not offset by the position
	local = surfaceLocal(s);
    if isfield(s,'position')
        position = s.position;
    else
        position = [0,0,0];
    end
    
    uv = (local * (xyz)')';
end