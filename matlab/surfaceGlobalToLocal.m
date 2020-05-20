function uv = surfaceGlobalToLocal(s,xyz)
	local = surfaceLocal(s);
    if isfield(s,'position')
        position = s.position;
    else
        position = [0,0,0];
    end
    
    uv = (local * (xyz - s.position)')';
end