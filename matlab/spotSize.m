function s=spotSize(rays)
%spotSize(rays) std dev of x and y positions

axes = rays.surface.local;
if size(axes,2) == 3
    ax1 = axes(1,:)';
    ax2 = axes(2,:)';
else
    ax1 = axes(:,1);
    ax2 = axes(:,2);
end

x = rays.position * ax1;
y = rays.position * ax2;
valid = rays.valid;
dx = std(x(valid));
dy = std(y(valid));
s = sqrt(dx * dy);
end

