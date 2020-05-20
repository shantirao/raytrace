function [ varargout ] = pareval( fn, varargin )
%pareval(fn, varargin) Evalute fn over the cell array. Result is a column cell array
%   Detailed explanation goes here
% if ~matlabpool('size')
%   matlabpool('open');
% end

N = max(cellfun(@(x)iscell(x)*numel(x),varargin));
L = length(varargin);
params = cell(N,L);
M = nargout;
result = cell(1,M);

for i=1:N, result{i} = cell(N,1); end;

for j=1:L
    if iscell(varargin{j})
        for i=1:N
            params{i,j} = varargin{j}{i};
        end
    else
        for i=1:N
            params{i,j} = varargin{j};
        end
    end
end


parfor i=1:N
    a = cell(1,M);
    [a{:}] = feval(fn,params{i,:});
    result{i}=a;
end

for j=1:M
    a = cell(N,1);
    for i=1:N
        a{i} = result{i}{j};
    end
    varargout{j} = a;
end

end

