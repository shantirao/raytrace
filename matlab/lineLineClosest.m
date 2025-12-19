function [points,dist] = lineLineClosest(p,d)
 % Solve for the line that is perpendicular to both
 % stable if the lines intersect intersection. returns NaN if lines are parallel
% https://en.wikipedia.org/wiki/Skew_lines
% p = [0 0 0; 0 1 0]; d = normr([1 1 1; 1 0 0]);
% p = [0 0 0; 0 1 0]; d = normr([1 0 0; 1 1 0]);
% p = [1 0 0; 0 1 1]; d = normr([0 1 0; 1 0 0]); %should be [1 1 0; 1 1 1]

 A = p(1,:);
 B = p(2,:);
 a = d(1,:);
 b = d(2,:);

 n = cross(a,b);
 if n*n' == 0
   points = NaN; return;
 endif

 n1 = cross(a,n);
 n2 = cross(b,n);

 dist = [(B-A)*n2' / (a*n2'); (A-B)*n1' / (b*n1')];
 points = p + dist.*d;
% n = cross(d(1,:), d(2,:));
% if n*n' == 0
%   points = NaN; return;
% endif
%
% n1 = cross(d(1,:),n);
% n2 = cross(d(2,:),n);
%
% dist = [(p(2,:)-p(1,:))*n2' / (d(1,:)*n2'); (p(1,:)-p(2,:))*n1' / (d(2,:)*n1')];
% point = p + dist.*d;

