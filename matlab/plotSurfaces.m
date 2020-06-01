function plotSurfaces(trace)
%plotSurfaces(start,finish)

xlabel('x');ylabel('y');zlabel('z');

if ~iscell(trace)
    trace = {trace};
end

hold on;
for i=2:numel(trace)
    if max(trace{i}.position,[],1) ~= [0,0,0]       
        v = trace{i}.valid;
        if isfield(trace{i},'surface') && isfield(trace{i},'local')
            x = trace{i}.position(v,:) * trace{i}.local(1,:)';
            y = trace{i}.position(v,:) * trace{i}.local(2,:)';
        else
            x = trace{i}.position(v,1);
            y = trace{i}.position(v,2);
        end
        if ~isempty(x) 
            tri = delaunay(x,y);
            z = trace{i}.position(v,3);
    %         v = v(tri(:,1)) & v(tri(:,2)) & v(tri(:,3));
    %         C = repmat([1 0 0],size(x)); %uint8(v(v));
            %delete unused areas

            h = trisurf(tri,x,y,z,'EdgeColor','none','FaceVertexCData',[1,0,0]); %,'CDataMapping','direct');
    %         colormap([1 0 0]);
            alpha(.5)
        end
    end
end
hold off;
end

