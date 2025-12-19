function [m,v] = nanmean (X, varargin)
  v = ~isnan(X);
  n = sum (v, varargin{:}); % varargin choose an axis for sum()
  n(n == 0) = NaN; %avoid divide by zero
  X(~v) = 0;
  m = sum (X, varargin{:}) ./ n;
end

