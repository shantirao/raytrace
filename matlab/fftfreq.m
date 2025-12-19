function f = fftfreq(n, d)
%  Given a window length n and a sample spacing d:
%
%f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
%f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
  if nargin < 2
    d = 1.0;
  end
  if mod(n,2) %odd
    f = [0, 1:(n-1)/2-1, -(n-1)/2:-1] / (d*n);
  else %even
    f = [0, 1:n/2-1, -n/2:-1] / (d*n);
  end

end
