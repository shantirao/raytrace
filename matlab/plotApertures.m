function [points, xy] = plotApertures(s,varargin)
% plotApertures(surface, is3D, {other plot specifications}
is3D = false;
if numel(varargin) && varargin{1} == true
    is3D = true;
end
if iscell(s)
    [points, xy] = map(@(x)plotApertures(x,varargin{:}),s);
%     points = s; %take advantage of copy-on-write to replace elements
%     for i=1:numel(s)
%         points{i} = plotApertures(s{i},varargin{:});
%     end
else
    style = {};
    points = {};
    xy = [];
    hold on;
    if isfield(s,'center') %offset to aperture geometry
        c = s.center(1:2); % back to the tangent plane projection
    else
        c = [0,0];
    end
    if is3D % aperture center (not necessarily origin for clipping)
        pos = surfaceLocalToGlobal(s,[0,0,elevation(s,[0,0])]);
    else
        pos = c;
    end
    if isfield(s,'style')
        style = s.style;
    end
    % calculate edge points x and y relative to the segment center
    if isfield(s,'diameter') && ~isfield(s,'aperture')
        s.aperture.type = 'circle';
        s.aperture.radius = s.diameter / 2;
    end
    if isfield(s,'aperture')
        ap = s.aperture;
        if isstruct(ap)
            if strcmp(ap.type,'pie')
                e = normr(ap.edges);
                if isfield(ap,'origin')
                  offset = ap.origin - c;
                else
                  offset = -c;
                end
% todo: start with a circle, then pull points in that exceed the edges or the annulus
%       project to edge to find the nearest point, and use the sign of x-product to determine if it's inside or outside.
Nsteps = 18;
                if isfield(ap,'gap')
                  d1 = ap.gap*normr(([0 1;-1 0]*e(1,:)')');%perpendicular to edges 1 and 2
                  d2 = ap.gap*normr(([0 -1;1 0]*e(2,:)')');
                else
                  d1 = d2 = 0;
                end
                angles = (1:Nsteps-1)*(acos(dot(e(1,:),e(2,:)))/Nsteps);
                xy = []; %ap.radius(1)*e(2,:)]+d2;
                for i=0:Nsteps % out on edg 2
                  xy(end+1,:)  = (ap.radius(1)*(Nsteps-i) + i*ap.radius(2))/Nsteps*e(2,:)+d2;
                end
                 last = xy(end,:);
                for a=angles  % todo: test if angle is outboard of the edge
                  v = ap.radius(2)*([cos(a) -sin(a);sin(a) cos(a)]*e(2,:)')'; %+d2
                  if numel(d2) > 1 && dot(d2,v - last) < 0, continue; end %must point away from the gap
                  xy(end+1,:) = v;
                end
%                  xy(end+1,:) = ap.radius(2)*e(1,:)+d1;
                for i=0:Nsteps % back on edge 1
                  xy(end+1,:) = (ap.radius(2)*(Nsteps-i) + i*ap.radius(1))/Nsteps*e(1,:)+d1;
                end
%                angles = fliplr(angles);
                last = xy(end,:);
                for a=angles % inner bit
                  v = ap.radius(1)*([cos(a) sin(a);-sin(a) cos(a)]*e(1,:)')';
                  if numel(d1) > 1 && dot(d1,v - last) < 0, continue; end
                  xy(end+1,:) = v;
                end
%                xy(end+1,:) = ap.radius(1)*e(2,:)+d2;
                if isfield(ap,'center') && isfield(ap,'parentRadius')
                  toofar = sum((xy - ap.center).^2,2) > ap.parentRadius.^2;
                  xy(toofar,:) = ap.center + ap.parentRadius * normr(xy(toofar,:)-ap.center);
                end
                xy = unique(xy,'rows','stable');
                xy(end+1,:) = xy(1,:); %close the shape
                x = xy(:,1) + offset(1);
                y = xy(:,2) + offset(2);
            elseif strcmp(ap.type,'circle')
                th = 0:pi/50:2*pi;
                x = ap.radius * cos(th');
                y = ap.radius * sin(th') ;
            elseif strcmp(ap.type,'annulus')
                th = 0:pi/50:2*pi;
                for j=1:numel(ap.radius)
                    x(:,j) = ap.radius(j) * cos(th');
                    y(:,j) = ap.radius(j) * sin(th') ;
                end
            elseif strcmp(ap.type,'rectangle')
                x1 = ap.bounds(1,1);
                x2 = ap.bounds(2,1);
                y1 = ap.bounds(1,2);
                y2 = ap.bounds(2,2);
                x = [x1 x1 x2 x2 x1]';
                y = [y1 y2 y2 y1 y1]';
            elseif strcmp(ap.type,'polygon')
                x = ap.bounds(:,1);
                y = ap.bounds(:,2);
                if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
                    x(end+1)=x(1);
                    y(end+1)=y(1);
                end
            end % ap.type

        else % isstruct -- legacy hack
          M = numel(ap);
          scale = max(max(ap(:)),max(ap(:))-min(ap(:)))/4;
          if M == 4 % rectangle [x1 y1 ; x2 y2]
              x1 = ap(1,1);
              x2 = ap(2,1);
              y1 = ap(1,2);
              y2 = ap(2,2);
              x = [x1 x1 x2 x2 x1]';
              y = [y1 y2 y2 y1 y1]';
              %plot(x + c(1),y+c(2));
          elseif M == 3 % ellipse
              error('Invalid aperture shape')
          elseif M < 3 % annulus or circle
              th = 0:pi/36:2*pi;
              for j=1:M
                  x(:,j) = ap(j) * cos(th');
                  y(:,j) = ap(j) * sin(th') ;
  %               %plot(x + c(1),y+c(2));
              end
          else % polygon
              x = ap(:,1);
              y = ap(:,2);
              if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
                  x(end+1)=x(1);
                  y(end+1)=y(1);
              end
%                 x = x;
%                 y = y;
  %             plot(x,y);
          end
        end
        xx = x;
        yy = y;
        for j=1:size(xx,2) %multiple shapes in the case of an annulus
            x = xx(:,j);
            y = yy(:,j);
            if (is3D)
                [z,~,xyz] = elevation(s,[x y]); %offsets center
%                xyz = surfaceLocalToGlobal(s,[x y z]); %offsets center
                points{end+1}=xyz;
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{2:end},style{:});
            else
                % surfaceLocalToTangent
                plot(x+c(1),y+c(2),varargin{2:end},style{:});
                points{end+1}=[x+c(1),y+c(2),zeros(size(y))];
            end
        end
        if isempty(xy)
          xy = [x(:) y(:)];
        end
    end
%     if (is3D)
%         if isfield(s,'local')
%             pos = surfaceLocalToGlobal(s,[0,0,elevation(s,[0,0])]);
%         else
%             xyz = s.position;
%         end
%     end
    if isfield(s,'referenceCenter')
      pos = s.referenceCenter;
    endif
    if ~isempty(points) && isfield(s,'reference') %reference axes are in global coordinates
        scale = max(range(points{1},1))/5;
        if (is3D)
            u = s.reference;
            plotVector(pos,scale*u(1,:),'r');
            plotVector(pos,scale*u(2,:),'g');
            plotVector(pos,scale*u(3,:),'b');
        else
            u = s.reference(:,1:2);
            plot(pos(1)+scale*[0 u(1,1)],pos(2)+scale*[0 u(1,2)],'r');
            plot(pos(1)+scale*[0 u(2,1)],pos(2)+scale*[0 u(2,2)],'g');
            plot(pos(1)+scale*[0 u(3,1)],pos(2)+scale*[0 u(3,2)],'b');
        end
    end
    if isfield(s,'name')
        if (is3D)
            text(pos(1),pos(2),pos(3),s.name,'HorizontalAlignment','center');
        else
            text(pos(1),pos(2),s.name,'HorizontalAlignment','center');
        end
    end
    hold off;
    if isfield(s,'segments') && numel(s.segments)
        % todo: map segment local coordinates onto parent local coordinates
        [p,sxy]=plotApertures(s.segments,varargin{:}); % passes is3D as the 2nd argument
        points = {points,p{:}};
        xy = {xy, sxy{:}};
    end
    if numel(points) == 1
      points = points{1};
    end
end
end
