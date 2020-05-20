function [ varargout ] = map(f, x )
%map calls f(x[i]) for each i in x, returns a cell array

L = numel(x);
N = nargout;
c = cell(size(x));
threshold = 6;

if iscell(x)
        for i=1:L
            if N>1
                t = cell(1,N);
                [t{:}]=f(x{i});
                c{i}=t;
            else
                c{i}=f(x{i});
            end
        end
else
        for i=1:L
            y = x(i);
            if N>1
                t = cell(1,N);
                [t{:}]=f(y);
                c{i}=t;
            else
                c{i}=f(y);
            end
        end        
end

if N > 1
    for i=1:N
        t = cell(size(x)); % output cell array
        for j=1:L
            y = c{j};
            t{j}=y{i};
        end
        varargout{i}=t;
    end
else
    varargout{1} = c;
end

end

