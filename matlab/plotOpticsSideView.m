function h=plotOpticsSideView(trace,axes,varargin)
%plotSideView(trace,[x-axis;y-axis],plotoptions)

if ~iscell(trace)
    trace = {trace};
end

if nargin < 2
    ax1 = [1;0;0];
    ax2 = [0;0;1];
elseif size(axes,2) == 3
    ax1 = axes(1,:)';
    ax2 = axes(2,:)';
else
    ax1 = axes(:,1);
    ax2 = axes(:,2);
end

x = cellfun(@(r)r.position * ax1 ,trace,'UniformOutput',false);
% x = map(@(r)r.position(:,1),raytrace);
x = horzcat(x{:});

y = cellfun(@(r)r.position * ax2,trace,'UniformOutput',false);
y = horzcat(y{:});

% valid = cellfun(@(r)r.valid,trace,'UniformOutput',false);
% valid = and(valid{:});
valid = trace{end}.valid;

ax3 = cross(ax1,ax2);
z = trace{end}.position * ax3;
z = (z == 0);
if sum(z) > 0
    valid = valid & z;
end

h=plot(x(valid,:),y(valid,:),varargin{:});

end

