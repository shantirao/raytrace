function local = surfaceLocal(s)
    if isstruct(s)
        if isfield(s,'local') && numel(s.local) >= 6
            local = s.local;
        else
            local = doSurfaceLocal(s.direction);
        end
    else
        local = doSurfaceLocal(s);
    end
end

function local = doSurfaceLocal(d)
    if abs(dot([0 1 0],d)) > 0.9
        local = normr(cross(d,[1 0 0]));
        local = [normr(cross(local,d));local];
    else
        local = normr(cross([0 1 0],d));
        local = [local; normr(cross(d,local))];
    end
end