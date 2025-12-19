function s=spotSize(rays)
%spotSize(rays) std dev of x and y positions

if ~isfield(rays.surface, 'local')
  avg = mean(rays.position(rays.valid,:),1);
  dist = rays.position(rays.valid,:) - avg;
  s = sqrt(mean(sum(dist.^2,2),1));
else
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
  dx = std(x(rays.valid));
  dy = std(y(rays.valid));
  s = sqrt(dx^2 + dy^2);
  end
end
