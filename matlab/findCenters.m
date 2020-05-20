function surfaces = findCenters(surfaces,source)
    options.aperture=false;
    options.segments=false;
    options.negative=true;
    [~,~,~,elev,trace] = raytrace(columnSource(source,0),surfaces,options);
    % trace the chief ray only
    N = numel(trace);
    for i=2:N
        if ~isfield(surfaces{i-1},'local')
            surfaces{i-1}.local = surfaceLocal(surfaces{i-1});
        end
        if ~isfield(surfaces{i-1},'center')
            surfaces{i-1}.center = surfaceTangentToLocal(surfaces{i-1},trace{i}.projection);
        end
        surfaces{i-1}.center(3) = elev{i};
    end
end
%         
% function local = makelocal(s)
% local = normr(cross([0 1 0],s.direction));
% local = [local; normr(cross(s.direction,local))];
% end