function rays = trimRays(rays) 
%removes invalid rays and 
if isfield(rays,'valid')
    v = rays.valid;
    fn = fieldnames(rays);
    for i=1:numel(fn)
        d = getfield(rays,fn{i});
        if size(d,1) == size(v,1)
            rays = setfield(rays,fn{i},d(v,:));
        end
    end    
    rays.N = numel(v(v));
end
end