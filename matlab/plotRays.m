function plotRays(trace,varargin)
%plotRays3({rays},options)

if ~iscell(trace)
    error('plotRays wants a cell array of ray positions')
end

x = cellfun(@(r)r.position(:,1),trace,'UniformOutput',false);
x = horzcat(x{:});

y = cellfun(@(r)r.position(:,2),trace,'UniformOutput',false);
y = horzcat(y{:});

z = cellfun(@(r)r.position(:,3),trace,'UniformOutput',false);
z = horzcat(z{:});

for i=1:numel(trace)
    if isfield(trace{i},'valid')
        x(~trace{i}.valid,i) = NaN;
        y(~trace{i}.valid,i) = NaN;
        z(~trace{i}.valid,i) = NaN;
    end
end
N = size(x,1);
downsample = ceil(N / 512);

sample=1:downsample:N;

h=plot3(x(sample,:)',y(sample,:)',z(sample,:)',varargin{:});

end

