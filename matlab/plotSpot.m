function [h,x,y]=plotSpot(rays,varargin)
%plotSpot(rays,plotoptions) plots spots in tangent plane coordinates
if iscell(rays)
    [M,N]=size(rays);i=1;
    for i=1:numel(rays)
        subplot(M,N,i);
        r = rays{i};
        plotSpot(r,varargin{:});
        if isfield(r,'surface') && isfield(r.surface,'name')
            title(r.surface.name);
        end
        i=i+1;
    end
else
if isfield(rays,'surface') %&& isfield(rays.surface,'local')
    xy = surfaceGlobalToLocal(rays.surface,rays.position);

    x = xy(:,1);
    y = xy(:,2);
%     axes = surfaceLocal(rays.surface);
%     origin = rays.surface.position;    
%     if isfield(rays.surface,'center')
%     end
else
    offset = [0,0];
    axes = [1 0 0;0 1 0];
    origin = [0 0 0];
    
    if size(axes,2) == 3
        ax1 = axes(1,:)';
        ax2 = axes(2,:)';
    else
        ax1 = axes(:,1);
        ax2 = axes(:,2);
    end

    pos = bsxfun(@minus,rays.position,origin);
    x = pos * ax1 - offset(1);
    y = pos * ax2 - offset(2);
end
colors = [];

% a = max(std(diff(x(valid))),std(diff(y(valid))));
a = max(x)-min(x);

if isfield(rays,'valid')
    valid = rays.valid;
    x = x(valid);
    y = y(valid);
    if isfield(rays,'segment')
        segment = rays.segment(valid);        
    end
else   
    if isfield(rays,'segment')
        segment = rays.segment;
    end
end
if isfield(rays,'segment')
    colors = varycolor(1+max(segment(:)));
    colors = colors(segment,:);
end
if numel(x) > 1024
    n = datasample(1:numel(x), 1024,'Replace',false);
    x = x(n);
    y = y(n);
    if size(colors,1)>numel(n),     colors = colors(n,:); end
end
if isempty(colors) || nargin>1
    h=scatter(x,y,'filled',varargin{:}); %a,varargin{:});
else
    h=scatter(x,y,16,colors,'filled'); %a,varargin{:});
end    
if isfield(rays,'units')
    xlabel(rays.units);
    ylabel(rays.units);
end
if isfield(rays,'surface') && isfield(rays.surface, 'aperture')
    plotApertures(rays.surface);
end
axis equal;

%% debug
% hold on;
% j =1;
% for i=1:rays.N
%     if rays.valid(i)
%         text(x(j),y(j),num2str(j));
%         j = j+1;
%     end
% end
% hold off;
end
end