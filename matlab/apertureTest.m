function [inAP, uv] = apertureTest(s, projection)
  % returns uv in the tangent plane projection of the local optic center

if isfield(s,'diameter') && ~isfield(s,'aperture')
  s.aperture = s.diameter / 2;
end
if isfield(s,'local')
  if isfield(s,'center')
    uv = projection * s.local(1:2,:)' - s.center(1:2);
  else
    uv = projection * s.local(1:2,:)';
  end
elseif isfield(s,'center') && numel(s.center) < 3 % 2d center but no local CS
  error('Define a local coordinate system if the aperture has a center point')
%    elseif isfield(s,'center') %center is in 3D
%      uv = projection - center;
else
  uv = projection;
end

if isstruct(s.aperture)
   if isfield(s.aperture,'origin')
     uv = uv - s.aperture.origin;
   end

   if strcmp(s.aperture.type,'pie') % gap, edges, origin, radius(inner,outer), parentRadius, parentCenter
      % [N x 3] * [3 x 2] - [1 x 2]
      % pie slice is with respect to a center and
      % winding order is u*y - v*x. [y1 y2; -x1 -x2] = [0 1;-1 0] * [x1 x1; y2 y2]'
      r2 = sum(uv.^2,2); % circular aperture; use winding order, same side of origin
      edges = normr(s.aperture.edges);
      winding = uv(:,1) * edges(:,2)' - uv(:,2) * edges(:,1)' + [-1 1].*s.aperture.gap(:)';% u y = v x
      side = sign(winding);
      q = sign(uv * edges) >= 0; % in the right quadrant?
      % opposite signs or on the line (zero) or at a vertex.
      inAP = ((side(:,1) ~= side(:,2))|(side(:,1)==0 & side(:,2)==0)) & q(:,1) & q(:,2) & s.aperture.radius(1)^2 <= r2 & r2 <= s.aperture.radius(2)^2; % inner, outer radius
%      inAP = ((side(:,1) ~= side(:,2))|(side(:,1)==0)) & dir(:,1) & dir(:,2) & ...
%          s.aperture.radius(1)^2 <= r2 & r2 <= s.aperture.radius(2)^2; % inner, outer radius
      if isfield(s.aperture,'parentRadius') && isfield(s.aperture,'center')
          %w.r.t. segment shape, not parent aperture
          uv = projection * s.local' - s.aperture.center(1:2);
          inAP = and(inAP, sum(uv.^2,2) < s.aperture.parentRadius.^2);
      end
  elseif strcmp(s.aperture.type,'circle')
      r2 = sum(uv.^2,2); % circular aperture
      inAP = r2 <= s.aperture.radius^2;
  elseif strcmp(s.aperture.type,'annulus')
      r2 = sum(uv.^2,2); % circular aperture
      inAP = s.aperture.radius(1)^2 <= r2 & r2 <= s.aperture.radius(2)^2;
  elseif strcmp(s.aperture.type,'rectangle')
      % s.aperture.bounds = |u1 u2|
      %                           |v1 v2|
      inAP = s.aperture.bounds(1) <= uv(:,1) & uv(:,1) <= s.aperture.bounds(3)& ...
          s.aperture.bounds(2) <= uv(:,2) & uv(:,2) <= s.aperture.bounds(4);
  elseif strcmp(s.aperture.type,'polygon')
      x = s.aperture.bounds(:,1);
      y = s.aperture.bounds(:,2);
      if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
          x(end+1)=x(1);
          y(end+1)=y(1);
      end
      [in,on] = inpolygon(uv(:,1),uv(:,2),x,y);
      inAP = in|on;
  else
      error(['invalid aperture type ' s.aperture.type]);
  end %aperture is a structure
elseif isfield(s,'aperture') && ~isempty(s.aperture) % center, local coordinate sys, radius|rectangle|ellipse|polygon
    M = numel(s.aperture);
    if M == 4 % rectangle
        % |u1 v1|  contains [u v]
        % |u2 v2|
        inAP = s.aperture(1) <= uv(:,1) & uv(:,1) <= s.aperture(2)& ...
            s.aperture(3) <= uv(:,2) & uv(:,2) <= s.aperture(4);
    elseif M == 3 % ellipse
        error('Invalid aperture shape')
    elseif M < 3
        r2 = sum(uv.^2,2); % circular aperture
        if M == 2 % [inner outer]
            inAP = s.aperture(1)^2 <= r2 & r2 <= s.aperture(2)^2;
        elseif M == 1 % outer
            inAP = r2 <= s.aperture^2;
        else
            error([label ' aperture must have 1, 2, or 4 elements']);
        end
    else % polygon
        x = s.aperture(:,1);
        y = s.aperture(:,2);
        if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
            x(end+1)=x(1);
            y(end+1)=y(1);
        end
        [in,on] = inpolygon(uv(:,1),uv(:,2),x,y);
        inAP = in|on;
    end
end
