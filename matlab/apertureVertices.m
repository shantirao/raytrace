function points = apertureVertices(s,n,varargin)
% Returns a list of vertices. 4 points for a pie slice, default 6 points for a circle
% use surfaceLocalToGlobal(s,apertureVertices(s)) to find 3D coordinates
if iscell(s)
    points = cell(numel(s),1);
    for i=1:length(s)
        points{i} = plotApertures(s{i},varargin{:});
    end
else
    hold on;
    if isfield(s,'center')
        c = s.center(1:2); % back to the tangent plane projection
    else
        c = [0,0]; 
    end
    
    points = {};
    scale = 100;
    if isfield(s,'aperture')
        ap = s.aperture;
        if isstruct(ap)
            if strcmp(ap.type,'pie')
                e = normr(ap.edges);
                offset = ap.origin - c; % parent origin, even though 
                e1 = e(1,:);
                e2 = e(2,:);
                n1 = [0 1;-1 0]*e1';
                n2 = [0 1;-1 0]*e2';
                points = [ap.radius(1)*e1-n1'*ap.gap(1);  ...
                    ap.radius(2)*e1-n1'*ap.gap(1); ...
                    ap.radius(2)*e2+n2'*ap.gap(2); ...
                    ap.radius(1)*e2+n2'*ap.gap(2)];
                points = bsxfun(@plus,points,offset);
            end                    
        else
            M = numel(ap);
            scale = max(max(ap(:)),max(ap(:))-min(ap(:)))/4;
            if M == 4 % rectangle
                x1 = ap(1,1);
                x2 = ap(2,1);
                y1 = ap(1,2);
                y2 = ap(2,2);
                x = [x1 x1 x2 x2 x1]';
                y = [y1 y2 y2 y1 y1]';
                points = [x y];
                %plot(x + c(1),y+c(2));
            elseif M == 3 % ellipse
                error('Invalid aperture shape')
            elseif M < 3 % annulus or circle
                if nargin < 2
                    n = 6;
                end
                th = (0:(n-1))/n * 2*pi;
                points = max(ap) *[ cos(th') sin(th')] ;
            else % polygon
                points = ap;
            end
        end
    else
        error('No aperture defined');
    end
%     points = bsxfun(@plus,points,c);
    z = elevation(s,points);
    points = [points z];
end
end