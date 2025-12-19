function in = pnpoly(x, y, xv, yv) % x,y are points, xv,yv are vertices
%https://wrfranklin.org/pnpoly
%slightly slower than the built-in inpolygon function
  nvert = numel(xv);
  in = false(size(x));
  c = in;
  j = nvert;
  for i=1:nvert
    a = ((yv(i)>y) != (yv(j)>y));
    b = a & (yv(j) > yv(i));
    c(b) =((x(b) - xv(i)) * (yv(j)-yv(i))) < (xv(j)-xv(i)) * (y(b) - yv(i));
    b = ~b;
    c(b) = ((x(b) - xv(i)) * (yv(i)-yv(j))) < (xv(j)-xv(i)) * (y(b) - yv(i));
    b = a&c;
    in(b) = ~in(b);

%  for i=1:nvert
%    if ((yv(i)>y) != (yv(j)>y))
%      if yv(j) > yv(i)
%        if ((x - xv(i)) * (yv(j)-yv(i))) < (xv(j)-xv(i)) * (y - yv(i))
%          in = ~in;
%        end
%      else % flip < to > because dividing by negative number. or swap i and j
%        if ((x - xv(i)) * (yv(i)-yv(j))) < (xv(j)-xv(i)) * (y - yv(i))
%          in = ~in;
%        end
%      end
%    end
    j = i; %next vertex in loop
  end
end
%% original relies on conditional evaluation to not divide by zero
% int pnpoly(int nvert, float *xv, float *yv, float testx, float testy)
%{
%  int i, j, c = 0;
%  for (i = 0, j = nvert-1; i < nvert; j = i++) {
%    if ( ((yv[i]>testy) != (yv[j]>testy)) &&
%	 (testx < (xv[j]-xv[i]) * (testy-yv[i]) / (yv[j]-yv[i]) + xv[i]) )
%       c = !c;
%  }
%  return c;
%}

%!test
%! pnpoly(.9,.9, [-1, -1, 1, 1], [-1, 1, 1, -1]) %true
%! pnpoly(0,0, [-1, -1, 1, 1], [-1, 1, 1, -1]) %true
%! pnpoly(2,0, [-1, -1, 1, 1], [-1, 1, 1, -1]) %false
%! pnpoly ([.95, 0, 2], [.95, 0, 0], [-1, -1, 1, 1], [-1, 1, 1, -1])
%! assert (in, [true, true, false]);
%! assert (on, [true, false, false]);


%!demo
%! xv = [ 0.05840, 0.48375, 0.69356, 1.47478, 1.32158, ...
%!        1.94545, 2.16477, 1.87639, 1.18218, 0.27615, ...
%!        0.05840 ];
%! yv = [ 0.60628, 0.04728, 0.50000, 0.50000, 0.02015, ...
%!        0.18161, 0.78850, 1.13589, 1.33781, 1.04650, ...
%!        0.60628 ];
%! xa = [0:0.1:2.3];
%! ya = [0:0.1:1.4];
%! [x,y] = meshgrid (xa, ya);
%! [in,on] = inpolygon (x, y, xv, yv);
%! inside = in & ! on;
%!
%! clf;
%! plot (xv, yv);
%! hold on;
%! plot (x(inside), y(inside), "og");
%! plot (x(! in), y(! in), "sm");
%! plot (x(on), y(on), "^b");
%! hold off;
%! disp ("Green circles are inside polygon, magenta squares are outside,");
%! disp ("and blue triangles are on the boundary.");

%!demo
%!  xv = [ 0.05840, 0.48375, 0.69356, 1.47478, 1.32158, ...
%!         1.94545, 2.16477, 1.87639, 1.18218, 0.27615, ...
%!         0.05840, 0.73295, 1.28913, 1.74221, 1.16023, ...
%!         0.73295, 0.05840 ];
%!  yv = [ 0.60628, 0.04728, 0.50000, 0.50000, 0.02015, ...
%!         0.18161, 0.78850, 1.13589, 1.33781, 1.04650, ...
%!         0.60628, 0.82096, 0.67155, 0.96114, 1.14833, ...
%!         0.82096, 0.60628];
%! xa = [0:0.1:2.3];
%! ya = [0:0.1:1.4];
%! [x, y] = meshgrid (xa, ya);
%! [in, on] = inpolygon (x, y, xv, yv);
%! inside = in & ! on;
%!
%! clf;
%! plot (xv, yv);
%! hold on;
%! plot (x(inside), y(inside), "og");
%! plot (x(! in), y(! in), "sm");
%! plot (x(on), y(on), "^b");
%! hold off;
%! disp ("Green circles are inside polygon, magenta squares are outside,");
%! disp ("and blue triangles are on the boundary.");

