function [surfaces,aperture] = makeRCTelescope(apertureDiameter,Fnumber,secondaryDistance,focusDistance)
%focus distance is location of the focal plane w.r.t. m1
%% Classic Cass design https://en.wikipedia.org/wiki/Cassegrain_reflector

% also http://www.telescope-optics.net/two-mirror.htm
focalLength = Fnumber * apertureDiameter;

% if where the pieces go is constrained and the m1 focal length can change
backFocalLength = secondaryDistance + focusDistance;
primaryFocalLength =  ( focalLength * secondaryDistance)/(-backFocalLength + focalLength - secondaryDistance);

% or maybe you have a particular f/number in mind for the  primary, 
% because of coating polarization concerns
% primaryFocalLength = Fnum1 * apertureDiameter; % focal length of the primary

% secondaryDistance = primaryFocalLength * (EFL - focusDistance)/(EFL + primaryFocalLength); 

Q = focalLength-backFocalLength-secondaryDistance; % distance from secondary vertex to focus

% m1 radius of curvature
R1 = 2*secondaryDistance*focalLength/(focalLength-backFocalLength); 
% m2 radius of curvature
R2 = 2*secondaryDistance*backFocalLength/Q; 

K1 = -1; % parabola
M = (focalLength-backFocalLength)/secondaryDistance; % secondary magnification
alpha = 1/2*sqrt(4*secondaryDistance*backFocalLength*M/((focalLength+backFocalLength*M-secondaryDistance*M)*(focalLength-backFocalLength-secondaryDistance)));
K2 = -1 - alpha - sqrt(alpha*(alpha+2));

%% useful for reporting; not needed in computing the prescription
mag = Q/(primaryFocalLength-secondaryDistance);
% EFL = primaryFocalLength * secondaryFocalLength / (primaryFocalLength - secondaryFocalLength - secondaryDistance)
Fnum1 = primaryFocalLength / apertureDiameter;

%% R-C https://en.wikipedia.org/wiki/Ritchey%E2%80%93Chr%C3%A9tien_telescope

  K1 = -1 - 2*backFocalLength/secondaryDistance/M^3;
  K2 = -1 - 2/(M-1)^3*(M*(2*M-1)+backFocalLength/secondaryDistance);

%% Prescription details

disp(sprintf('optic\tD   \tZ   \tR   \tK   \tEFL \tFnum'));
disp(sprintf('M1   \t0   \t0   \t%.5g\t%.5g\t%.5g\t%.5g',R1,K1,primaryFocalLength,Fnum1));
disp(sprintf('M2   \t%.5g\t%.5g\t%.5g\t%.5g\t%.5g',secondaryDistance,secondaryDistance,R2,K2));
disp(sprintf('PF   \t%.5g\t%.5g\t-\t-\t%.5g\t%.5g',backFocalLength,-focusDistance,focalLength,Fnumber));
disp(sprintf('M2 Magnification: %.5g',mag));
disp(sprintf('Focal length: %.5g',focalLength));

%%

aperture.position=  [0,0,secondaryDistance];
aperture.direction= [0,0,-1];
aperture.local = [-1 0 0; 0 1 0];
aperture.aperture = apertureDiameter/2;

m1.name = 'm1';
m1.type = 1; %reflective
m1.position =  [0,0,0];
m1.direction=[0,0,1];
m1.local = [1,0,0;0,1,0];
m1.aperture = apertureDiameter/2;
m1.cuy = 1 / R1;
m1.K = K1;
m1.n = 1;

m2.name = 'm2';
m2.type = 1; %reflective
m2.position = [0,0,secondaryDistance];
m2.direction=-m1.direction;
m2.local = [-1,0,0;0,1,0];
m2.cuy = -1 / R2;
m2.K = K2;
m2.n = 1;

pf.name = 'focus';
pf.type = 0;
pf.position = [0,0,-focusDistance];
pf.direction = m1.direction;
pf.local = m1.local;

surfaces ={m1,m2,pf};
end