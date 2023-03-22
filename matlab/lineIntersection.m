function [location,distance] = lineIntersection(position,direction)
%% Positions A,B, directions a,b, Solve (A-B + ax - by = 0)
% distance = [x; y], intersection = A+ax = B+by
%position = gather(position);
%direction = gather(direction);
offset = diff(position,1)';
A = [1 0; 0 -1]* direction;
distance = linsolve(A',offset);

location = position + diag(distance)*direction;
% faster than bsxfun(@times,distance,direction);

% if max(abs(diff(location,1))) / max(abs(location(:))) > 1e-7
%     error('intersection not within numerical tolerance')
% end
location = location(1,:);
end

%% test cases
% [p,d] = lineIntersection([0 0;0 4],[1 1;1 -2])
%  [p,d] = lineIntersection([0 0;0 4],[2 1;2 -3])
