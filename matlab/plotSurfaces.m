function plotSurfaces(trace)
%plotSurfaces(start,finish)

xlabel('x');ylabel('y');zlabel('z');

if ~iscell(trace)
    trace = {null,trace};
end

hold on;
for i=2:numel(trace)
    if max(trace{i}.position,[],1) ~= [0,0,0]       
        v = trace{i}.valid;
%         if isfield(trace{i},'surface') && isfield(trace{i},'local')
%             x = trace{i}.position(v,:) * trace{i}.local(1,:)';
%             y = trace{i}.position(v,:) * trace{i}.local(2,:)';
%             z = trace{i}.position(v,:) * trace{i}.surface.direction';
%         else
            x = trace{i}.position(v,1);
            y = trace{i}.position(v,2);
            z = trace{i}.position(v,3);
%         end
        if ~isempty(x) 
            if isfield(trace{i},'segment')
                s = trace{i}.segment(v);
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
        alpha(.5)
    end
end
