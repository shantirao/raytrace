function [ result ] = mcrun( N, fcn, varargin)
%mcrun Runs fcn(args{:}) N times, returns result as a cell array
result = cell(1,N);
if numel(varargin) == 1 && iscell(varargin{1})
    for i=1:N
        result{i} = feval(fcn,varargin{1}{:});
    end
else
    for i=1:N
    	result{i} = feval(fcn,varargin{:});
    end
end
% for i=1:N
%     result{i} = feval(fcn,args{:});
% end
end

