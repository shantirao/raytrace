function [ h ] = plotPolyLine(points,varargin)
%plotLine(start,finish,...)
points = [points;points(1,:)];
if size(points,2) == 3
    h = plot3(points(:,1),points(:,2),points(:,3),varargin{:});
else
    h = plot(points(:,1),points(:,2),varargin{:});
end

