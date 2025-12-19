function plotSurfaces(lightpath)
%plotSurfaces(start,finish)

xlabel('x');ylabel('y');zlabel('z');

if ~iscell(lightpath)
    lightpath = {null,lightpath};
end

hold on;
for i=2:numel(lightpath)
    if max(lightpath{i}.position,[],1) ~= [0,0,0]
        v = lightpath{i}.valid;
%         if isfield(lightpath{i},'surface') && isfield(lightpath{i},'local')
%             x = lightpath{i}.position(v,:) * lightpath{i}.local(1,:)';
%             y = lightpath{i}.position(v,:) * lightpath{i}.local(2,:)';
%             z = lightpath{i}.position(v,:) * lightpath{i}.surface.direction';
%         else
            x = lightpath{i}.position(v,1);
            y = lightpath{i}.position(v,2);
            z = lightpath{i}.position(v,3);
%         end
        if ~isempty(x)
            if isfield(lightpath{i},'segment')
                s = lightpath{i}.segment(v);
                mn = min(s(:));
                mx = max(s(:));
                for j=mn:mx
                    drawSurface(x,y,z,s(:) == j);
                end
            else
                drawSurface(x,y,z);
%                 tri = delaunay(x,y);
%                 h = trisurf(tri,x,y,z,'EdgeColor','none','FaceVertexCData',[1,0,0]); %,'CDataMapping','direct');
%                 alpha(.5)
            end
        end
    end
end
hold off;
end

function t = isCollinear(x)
 t = (max(x(:)) - min(x(:)) ) < 1e-14;
end

function drawSurface(x,y,z, select)
    if nargin > 3
        x=x(select);
        y=y(select);
        z=z(select);
    end
    if ~isCollinear(x) && ~isCollinear(y)
        tri = delaunay(x,y);
        h = trisurf(tri,x,y,z,'EdgeColor','none','FaceVertexCData',[1,0,0]); %,'CDataMapping','direct');
%        alpha(.5)
set(h,'edgecolor','none')
set(h,'facealpha',0.5)

    end
end
