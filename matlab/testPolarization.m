% polarization test

mrot = @(phi)[1 0 0 0;0 cos(2*phi) sin(2*phi) 0;0 -sin(2*phi) cos(2*phi) 0; 0 0 0 1];
horiz = [1 1 0 0]';
vert = [1 -1 0 0]';
diag1 = [1 0 1 0]';
diag2 = [1 0 -1 0]';

mrot(0) % should be identity
mrot(pi/2)*diag1 %should be diag2
mrot(-pi/4)*horiz %should be diag1

%%
% reflection of a circularly polarized beam off a corner cube should result
% in a flip in circularization

% rooftop prism in the XY plane
s1 = struct('position',[0 0 0], 'direction',normr([1,-1,0]),'type',1);
s2 = struct('position',[0 0 0], 'direction',normr([-1,-1,0]),'type',1);
s1.NK = 1+15i;
s2.NK = 1+15i;
% ray moves in the +Y direction, TE in the X direction, circularly
source = struct('position',[-1,-10,0],'direction',[0 1 0],'local',[1 0 0]);

% right hand circularly polarized
source.polarization=[1 0 0 1]; 
ray = raytrace(source,{s1,s2});
disp(ray.polarization);

% source.polarization=[1 0 0 1]; % vertically polarized
% ray = raytrace(source,{s1,s2});
% disp(ray.polarization);


%% fresnel matrix coefficients
% compare with http://people.physics.tamu.edu/mgao/calculator/my.jqplot.fresnel.html
n1=1; n2=1+15i; a=0; [rp, rs, m11, m12, m33, m34] =  fresnel(n1*cosd(10),n1,n2);
disp([rp, rs])
disp('   0.9821 + 0.1330i  -0.9829 - 0.1290i')
disp([m11, m12, m33, m34]);
disp('    0.9825 -0.00027 -0.9825 -0.0040')

%% multiple values?
[rp, rs, m11, m12, m33, m34] =  fresnel(n1*cosd([1:90]'),n1,n2);
