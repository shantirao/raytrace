function Geometry = telescopeGeometry(label,prescription,aperture,display)
%% find exit pupil and insert it before or after the focus 1in the prescription
% the last element of the prescription has to be a focus point.
% trace to the pupil relay instead of to the image
[pupilRelay,surfaces] = findExitPupil(prescription,aperture,display);
options.negative=true;  
% allow negative tracing from the focus back to a virtual pupil

if display
    savepng([label ' find pupil relay']);
end

% for calculating PSFs
% nBands = 4;
% bandWidth = 0.000080;
% bandpass = bsxfun(@plus,[0.000632,1],bandWidth*(-(nBands-1)/2:(nBands-1)/2)');
% bandpass(:,2) = 1/size(bandpass,1);


%% save geometry and prescription
Geometry = struct();
Geometry.prescription = surfaces;
Geometry.aperture = aperture;
if isfield(surfaces{1},'segments')
    Geometry.N = numel(surfaces{1}.segments);
else
    Geometry.N = 0;
end

%% calculate the extent of the secondary mirror
source = sourceColumn(aperture,1);
opt.segments= false;
opt.aperture=false;
opt.negative=true;
trace = raytrace(source,{surfaces{1:2}},opt);
[x,y] = spotPosition(trace{end});
surfaces{2}.center = [x,y];
[x,y] = spotExtent(trace{end});
surfaces{2}.aperture = [max(abs(x)),max(abs(y))]; 
%ellipse. Make it off-axis by moving the aperture

%%
Geometry.m1BeamLaunchers = cellfun(@(s)surfaceLocalToGlobal(s,apertureVertices(s)),surfaces{1}.segments,'UniformOutput',false);
Geometry.m2CornerCubes = surfaceLocalToGlobal(surfaces{2},apertureVertices(surfaces{2},12));

%% Display geometry
if display
    figure(2); clf; 
    subplot(1,3,1);plotApertures(surfaces{1},true);axis equal;
    title('M1 apertures')
    subplot(1,3,2);plotApertures(surfaces{2},true);axis equal; %segmented m2?
    plotDatums(Geometry.m2CornerCubes);
    title('M2 apertures and corner cube positions')
    subplot(1,3,3);
    plotApertures(surfaces{1},true);
    plotDatums(Geometry.m1BeamLaunchers);
    plotApertures(surfaces{2},true);
    plotDatums(Geometry.m2CornerCubes);

    axis equal;xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]');grid on;
    title({'Beam launcher and corner cube positions';'X = red, Y = green, Z = blue'});
    if numel(label)
    	savepng([label ' segment coordinate system']);
    end
end