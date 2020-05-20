function [ h ] = plotVector(origin,direction,varargin)
%plotVector(origin,direction,...)
if numel(origin) == 3
    h = plot3(origin(1) +[0;direction(:,1)],origin(2) +[0;direction(:,2)],origin(3) +[0;direction(:,3)],varargin{:});
else
    h = plot(origin(1) +[0;direction(:,1)],origin(2) +[0;direction(:,2)],varargin{:});

end

