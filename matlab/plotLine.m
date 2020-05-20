function [ h ] = plotLine(start,finish,varargin)
%plotLine(start,finish,...)
if numel(start) == 3
    h = plot3([start(1);finish(1)],[start(2);finish(2)],[start(3);finish(3)],varargin{:});
else
    h = plot([start(1);finish(1)],[start(2);finish(2)],varargin{:});
end

