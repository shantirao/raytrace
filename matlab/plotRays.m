function plotRays(trace,varargin)
%plotRays3(start,finish)

if ~iscell(trace)
    trace = {trace};
end

x = cellfun(@(r)r.position(:,1),trace,'UniformOutput',false);
% x = map(@(r)r.position(:,1),raytrace);
x = horzcat(x{:});

y = cellfun(@(r)r.position(:,2),trace,'UniformOutput',false);
y = horzcat(y{:});

z = cellfun(@(r)r.position(:,3),trace,'UniformOutput',false);
z = horzcat(z{:});

N = trace{1}.N;
downsample = ceil(N / 512);

sample=1:downsample:N;

h=plot3(x(sample,:)',y(sample,:)',z(sample,:)',varargin{:});

end

