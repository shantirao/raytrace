function f = plotDatums(p)
if iscell(p)
    f = {};
    for i=1:numel(p)
        f{i} = plotDatums(p{i});
    end
else
    hold on;
        f = scatter3(p(:,1),p(:,2),p(:,3),'ob');
        for j=1:size(p,1)
            text(p(j,1),p(j,2),p(j,3),num2str(j));
        end
    hold off;
end