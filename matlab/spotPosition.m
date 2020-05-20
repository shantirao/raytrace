function [x,y]=spotPosition(rays)
%spotSize(rays) std dev of x and y positions

axes = rays.surface.local;
if size(axes,2) == 3
    ax1 = axes(1,:)';
    ax2 = axes(2,:)';
else
    ax1 = axes(:,1);
    ax2 = axes(:,2);
end

x = mean(rays.position(rays.valid,:) * ax1);
y = mean(rays.position(rays.valid,:) * ax2);
end

