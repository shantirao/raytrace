function out = resizeAndCenter(in, rows, cols);
% pads or crops the center of the matrix into a new shape
   if nargin < 3
      cols = rows;
   end

   out = zeros(rows,cols);

   [m, n]  = size(in);

   height = min(m,rows);
   width = min(n,cols);

   i = max(0,floor(m/2-height/2));
   j = max(0,floor(n/2-width/2));

   y = min(height,m);
   x = min(width,n);
   u = max(0,floor(height/2-m/2));
   v = max(0,floor(width/2-m/2));

   out(i+1:i+y,j+1:j+x) = in(u+1:u+y,v+1:v+x);
 end

