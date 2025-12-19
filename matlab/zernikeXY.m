function [z, dxdy] = zernikeXY(c, xy, R)
% Nobody seems to agree, or know, which zernike term is which. So we'll just explicity name them and save it in a struct
% convert your freeform deformations to XY polynomials up to 8th order. Specify aperture.center if you want to offset to a specific location
% or just use an aperture deformation if you want the vertex to be the center
% xy should be normalized to a unit circle
x = xy(:,1) / R;
y = xy(:,2) / R;
N = size(x,1);
z = zeros(N,25);
dx = zeros(N,25);
dy = zeros(N,25);
i = 1;

%  SPS XYP 450.0
%  SCO NRADIUS 450.0; SCC NRADIUS 100
%Y   SCO C3 -1.77600216617074; SCC C3 100
%X2  SCO C4 17.07701268094119; SCC C4 100
%Y2  SCO C6 -16.90145809137451; SCC C6 100
%X2Y SCO C8 -2.61872493837637; SCC C8 100
%Y3  SCO C10 4.042909791206121; SCC C10 100
%X4  SCO C11 0.4965277804889697; SCC C11 100
%X2Y2  SCO C13 2.094638473092603; SCC C13 100
%Y4  SCO C15 -1.363343648699021; SCC C15 100
%X4Y  SCO C17 -0.1144979326662042; SCC C17 100
%X2Y3  SCO C19 -0.6027410401700329; SCC C19 100
%Y5  SCO C21 0.6330581321426945; SCC C21 100
%X6  SCO C22 0.4508135287733071; SCC C22 100
%X4Y2  SCO C24 0.6456644133048351; SCC C24 100
%X2Y4  SCO C26 0.3415378864774281; SCC C26 100
%Y6  SCO C28 -0.3432971024765515; SCC C28 100
%X6Y  SCO C30 -0.1478293724858129; SCC C30 100
%X4Y3  SCO C32 -0.1197216325716414; SCC C32 100
%X2Y5  SCO C34 0.2040448523141557; SCC C34 100
%Y7  SCO C36 0.1759371123999843; SCC C36 100
%X8  SCO C37 -0.224402981705639; SCC C37 100
%X6Y2  SCO C39 -0.7246066574109141; SCC C39 100
%X4Y4  SCO C41 -0.8274020819989081; SCC C41 100
%X2Y6  SCO C43 -0.37859611858763; SCC C43 100
%Y8  SCO C45 -0.05139771229399696; SCC C45 100

x2 = x.^2;
y2 = y.^2;
x3 = x2 .* x;
y3 = y2 .* y;
x4 = x2.^2;
y4 = y2.^2;
x5 = x2 .* x3;
y5 = y2 .* y3;
x6 = x3.^2;
y6 = y3.^2;
x7 = x3 .* x4;
y7 = y3 .* y4;
x8 = x4.^2;
y8 = y4.^2;

% compute up to 8th order
if isfield(c,'X')
  dx(:,i) = c.X;
  z(:,i++) = c.X * x;
end
if isfield(c,'Y')
  dy(:,i) = c.Y;
  z(:,i++) = c.Y * y;
end
if isfield(c,'X2')
  dx(:,i) = 2 * c.X2 *x;
  z(:,i++) = c.X2 * x2;
end
if isfield(c,'XY')
  dx(:,i) = c.XY * y;
  dy(:,i) = c.XY * x;
  z(:,i++) = c.XY * x .* y;
end
if isfield(c,'Y2')
dy(:,i) = 2*c.Y2 * y;
z(:,i++) = c.Y2 * y2;
end
if isfield(c,'X2Y')
dx(:,i) = 2*c.X2Y *x .*y;
dy(:,i) = c.X2Y * x2;
z(:,i++) = c.X2Y * x2 .* y;
end
if isfield(c,'X3')
dx(:,i) = 3 * c.X3 * x2;
z(:,i++) = c.X3 * x3;
end
if isfield(c,'Y3')
dy(:,i) = 3 * c.Y3 * y2;
z(:,i++) = c.Y3 * y3;
end
if isfield(c,'X4')
dx(:,i) = 4 * c.X4 * x3;
z(:,i++) = c.X4 * x4;
end
if isfield(c,'X2Y2')
dx(:,i) = 2 * c.X2Y2 * x .* y2;
dy(:,i) = 2 * c.X2Y2 * x2 .* y;
z(:,i++) = c.X2Y2 * x2 .* y2;
end
if isfield(c,'Y4')
dy(:,i) = 4 * c.Y4 * y3;
z(:,i++) = c.Y4 * y4;
end
if isfield(c,'X4Y')
dx(:,i) = c.X4Y * 4 * x3 .* y;
dy(:,i) = c.X4Y * x4;
z(:,i++) = c.X4Y * x4 .* y;
end
if isfield(c,'X2Y3')
dx(:,i) = c.X2Y3 * 2 * x .* y3;
dy(:,i) = c.X2Y3 * 3 * x2 .* y2;
z(:,i++) = c.X2Y3 * x2 .* y3;
end
if isfield(c,'X5')
dx(:,i) = c.X5 * 5 * x4  ;
z(:,i++) = c.X5 * x5 ;
end
if isfield(c,'Y5')
dy(:,i) = c.Y5 * 5 * y4  ;
z(:,i++) = c.Y5 * y5 ;
end
if isfield(c,'X6')
dx(:,i) = c.X6 * 6*x5 ;
z(:,i++) = c.X6 * x6 ;
end
if isfield(c,'X4Y2')
dx(:,i) = c.X4Y2 * 4 * x3 .* y2 ;
dy(:,i) = c.X4Y2 * 2 * x4 .* y  ;
z(:,i++) = c.X4Y2 * x4 .* y2 ;
end
if isfield(c,'X2Y4')
dx(:,i) = c.X2Y4 * 2 * x .* y4 ;
dy(:,i) = c.X2Y4 * 4 * x2 .* y3  ;
z(:,i++) = c.X2Y4 * x2 .* y4 ;
end
if isfield(c,'Y6')
dy(:,i) = c.Y6 * 6 * y5;
z(:,i++) = c.Y6 * y6 ;
end
if isfield(c,'X6Y')
dx += c.X6Y * 6 * x5 .* y ;
dy += c.X6Y * x6  ;
z(:,i++) = c.X6Y * x6 .* y ;
end
if isfield(c,'X4Y3')
dx(:,i) = c.X4Y3 * 4 * x3 .* y3  ;
dy(:,i) = c.X4Y3 * 3 * x4 .* y2  ;
z(:,i++) = c.X4Y3 * x4 .* y3 ;
end
if isfield(c,'X2Y5')
dx(:,i) = c.X2Y5 * 2 * x .* y5 ;
dy(:,i) = c.X2Y5 * 5 * x2 .* y4  ;
z(:,i++) = c.X2Y5 * x2 .* y5 ;
end
if isfield(c,'Y7')
dy(:,i) = c.Y7 * 7 * y6  ;
z(:,i++) = c.Y7 * y7 ;
end
if isfield(c,'X7')
dx(:,i) = c.X7 * 7 * x6  ;
z(:,i++) = c.X7 * x7 ;
end
if isfield(c,'X8')
dx(:,i) = c.X8 * 8 * x7 ;
z(:,i++) = c.X8 * x8 ;
end
if isfield(c,'X6Y2')
dx(:,i) = c.X6Y2 * 6 * x5 .* y2 ;
dy(:,i) = c.X6Y2 * 2 * x6 .* y  ;
z(:,i++) = c.X6Y2 * x6 .* y2 ;
end
if isfield(c,'X4Y4')
dx(:,i) = c.X4Y4 * 4 * x3 .* y4 ;
dy(:,i) = c.X4Y4 * 4 * x4 .* y3 ;
z(:,i++) = c.X4Y4 * x4 .* y4 ;
end
if isfield(c,'X2Y6')
dx(:,i) = c.X2Y6 * 2 * x .* y6 ;
dy(:,i) = c.X2Y6 * 6 * x2 .* y5 ;
z(:,i++) = c.X2Y6 * x2 .* y6 ;
end
if isfield(c,'Y8')
dy(:,i) = c.Y8 * 8 * y7 ;
z(:,i++) = c.Y8 * y8 ;
end

% numerical error from high dynamic range?
%symbolic version, very very slow for 100 digits of precision
%pkg load symbolic;
%sdx = double(sum(vpa(dx,100),2));
%sdy = double(sum(vpa(dy,100),2));
%sz = double(sum(vpa(z,100),2));

% sort to give the smaller numbers a chance to contribute
dx = sum(sort(dx,2),2);
dy = sum(sort(dy,2),2);
z = sum(sort(z,2),2);

dxdy = [dx dy] / R; % ( x = X / R, so  dx/dX = 1/R )

