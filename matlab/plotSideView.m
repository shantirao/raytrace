function h=plotSideView(trace,axes,varargin)
%plotSideView(trace,[x-axis;y-axis],plotoptions)

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
if ~iscell(trace)
    trace = {trace};
end

% x = cellfun(@(r)r.position * ax1 ,trace,'UniformOutput',false);
% x = horzcat(x{:});
% 
% y = cellfun(@(r)r.position * ax2,trace,'UniformOutput',false);
% y = horzcat(y{:});
% 
% valid = trace{end}.valid;

x = []; y=[];
for i=1:numel(trace)
   r = trace{i};
   x = [x, r.position * ax1];
   y = [y, r.position * ax2];
   if isfield(r,'valid')
       valid = r.valid;
       x(~valid,end) = NaN;
       y(~valid,end) = NaN;
   end
end

% ax3 = cross(ax1,ax2);
% z = trace{end}.position * ax3;
% z = (abs(z) <= 100);
% if sum(z) > 0
%     valid = valid & z;
% end

% h=plot(x(valid,:)',y(valid,:)',varargin{:});
h=plot(x',y',varargin{:});
axis equal;
end

