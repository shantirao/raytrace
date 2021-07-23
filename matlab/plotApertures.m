function points = plotApertures(s,varargin)

is3D = false;
if numel(varargin) && varargin{1} == true
    is3D = true;
end
if iscell(s)
    points = map(@(x)plotApertures(x,varargin{:}),s);
%     points = s; %take advantage of copy-on-write to replace elements
%     for i=1:numel(s)
%         points{i} = plotApertures(s{i},varargin{:});
%     end
else
    points = {};
    style = {};
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
                offset = ap.origin - c; % parent origin, even though 
%                 xy = [ap.radius(1)*e(1,:) ; ...
%                     ap.radius(2)*e(1,:) ; ...
%                     ap.radius(2)*e(2,:) ; ...
%                     ap.radius(1)*e(2,:) ; ...
%                     ap.radius(1)*e(1,:)];
                angles = 0:0.05:acos(dot(e(1,:),e(2,:)));
                xy = [ap.radius(1)*e(1,:) ;  ap.radius(2)*e(1,:)];
                %not computationally efficient but who cares?
                for a=angles
                    xy = [xy; ap.radius(2)*([cos(a) -sin(a);sin(a) cos(a)]*e(1,:)')'];
                end
                xy = [xy; ap.radius(2)*e(2,:) ;  ap.radius(1)*e(2,:)];
                angles = flip(angles);
                for a=angles
                    xy = [xy; ap.radius(1)*([cos(a) -sin(a);sin(a) cos(a)]*e(1,:)')'];
                end
                xy = [xy;ap.radius(1)*e(1,:)];
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
            else %legacy hack
            M = numel(ap);
            scale = max(max(ap(:)),max(ap(:))-min(ap(:)))/4;
            if M == 4 % rectangle
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
                th = 0:pi/50:2*pi;
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
                z = elevation(s,[x y]); %offsets center
                xyz = surfaceLocalToGlobal(s,[x y z]); %offsets center
                points{end+1}=xyz;
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),varargin{2:end},style{:});
            else
                % surfaceLocalToTangent
                plot(x+c(1),y+c(2),varargin{:},style{:});
                points{end+1}=[x+c(1),y+c(2),zeros(size(y))];
            end
        end
    end
%     if (is3D)
%         if isfield(s,'local')
%             pos = surfaceLocalToGlobal(s,[0,0,elevation(s,[0,0])]);
%         else
%             xyz = s.position;
%         end
%     end
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
        p=plotApertures(s.segments,varargin{:});
        points = {points{:},p{:}};
    end        
end
end