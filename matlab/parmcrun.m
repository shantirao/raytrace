function [ result ] = parmcrun( N, fcn, varargin)
%mcrun Runs fcn(args{:}) N times, returns result as a cell array
if ~matlabpool('size')
  matlabpool('open');
end

result = cell(1,N);
if numel(varargin) == 1 && iscell(varargin{1})
    parfor i=1:N
        result{i} = feval(fcn,varargin{1}{:});
    end
else
    parfor i=1:N
    	result{i} = feval(fcn,varargin{:});
    end
end
end

