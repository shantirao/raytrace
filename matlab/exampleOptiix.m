if exist('OCTAVE_VERSION', 'builtin') ~= 0
  pkg load statistics;
end
%% Design parameters. Dimensions in mm.
apertureDiameter = 3800;
Fnum = 6;
numRings = 2;
secondaryDistance = 4000;
backFocalLength = 2000; % distance behind m1 to prime focus, [mm]
segmentSize = 1500;
segmentSpacing = 50;
% mag = 29.17, not used
[prescription,aperture] = makeRCTelescope(apertureDiameter,Fnum,secondaryDistance,backFocalLength);
prescription{1}.segments = makeHexSegments(prescription{1},[0,0],segmentSize,segmentSpacing,pi/6,numRings);
prescription{1}.segments = {prescription{1}.segments{2:end}}; % no center
aperture.units = 'mm';
prescription{end}.display = 'nm';
label = 'OpTIIXRC';
%%
Geometry = telescopeGeometry(label,prescription,aperture,true);
save([label 'Geometry'],'Geometry');

%% print prescription
surfaceCenters={};
disp(sprintf('name\tcenter\tposition\tdirection\taperture\tR\tK'));
for i=1:numel(Geometry.prescription)
    p = Geometry.prescription{i};
    if isfield(p,'aperture'), ap = sprintf('%d,',p.aperture); else ap='';end
    if isfield(p,'cuy'), cuy = sprintf('%d,',p.cuy); else cuy='';end
    if isfield(p,'k'), k = sprintf('%d,',p.K); else k='';end
    c = surfaceLocalToGlobal(p,[0,0,0]);
    surfaceCenters{end+1}=c;
    disp(sprintf('%s\t[%g,%g,%g]\t[%g,%g,%g]\t[%g,%g,%g]\t[%s]\t%s\t%s',...
        p.name,c,p.position,p.direction,ap,cuy,k));
    if isfield(p,'segments')
        for j=1:numel(p.segments)
            s = p.segments{j};
            c = surfaceLocalToGlobal(s,[0,0,0]);
            disp(sprintf('%s\t[%g,%g,%g]',s.name,c));
        end
    end
end
%save Geometry;
save([label 'surfaceCenters'],'surfaceCenters');
%%
disp('M1 beam launchers')
for i=1:numel(Geometry.prescription{1}.segments)
    disp(Geometry.prescription{1}.segments{i}.name)
    disp(sprintf('[%g, %g, %g]\n',Geometry.m1BeamLaunchers{i}'))
end
disp('M2 Corner cubes')
    disp(sprintf('[%g, %g, %g]\n',Geometry.m2CornerCubes'))

%% test one degree of freedom
source = sourceColumn(Geometry.aperture,99,1);
trace = raytrace(source,prescription,options);
trace = trace(2:end); %skip the source rays
nomRays = trace{end};
[nomOPL,nomMask] = pupilOPL(nomRays,false);
[rmswfe,piston] = rmsWFE(nomOPL,nomMask);
perturbationScale = 1e-3*[1,1,1,.001,.001,.001];
perturbations = diag(perturbationScale); % 1 nm and 1 nanoradians
labels = {'Local \Deltax','Local \Deltay','Local \Deltaz','Local \thetaX','Local \thetaY','Local \thetaZ'};
unitLabels = {'nm/um','nm/um','nm/um','nm/urad','nm/urad','nm/urad'};
options.negative=true; %trace backward at end from focal plane to pupil relay surface for beamwalk reasons
scale = 1000;
labelScale = 1e6; %mm to nm
figure(40);clf;
for i=1:6
  subplot(2,3,i)
  sp = prescription;
  for j=1:6
  sp{1}.segments{j} = perturb(prescription{1}.segments{j},perturbations(i,:),true);
  end
  prays=lastCell(raytrace(source,sp,options));
  [opl,m] = pupilOPL(prays,false);

  pwf = opl - nomOPL;
  pwf(~nomMask(:))=NaN;
  imagesc(scale*pwf); axis image; colorbar;
  title(sprintf("%s %.3d %s",labels{i},labelScale*sqrt(mean(pwf(nomMask(:)).^2)),unitLabels{i}));
end
%title(sprintf('All perturbations %s rms/um',surfaces{1}.segments{j}.name,units));


%%
nPoints = 1024;
tic
[Sensitivity,surfaces] = toleranceAnalysis(label,Geometry,nPoints,true);
toc
%%
save([label 'Sensitivity'],'Sensitivity');

%%
figure(6);clf;
offsets = [0 0 2 0 0 0];

scale = 1/Sensitivity.perturbations(1);
for i=1:6
    subplot(2,3,i);
    wfe = Sensitivity.wfe{i,1};
    dz = offsets(i);
    for j=2:Sensitivity.N
        wfe = wfe + Sensitivity.wfe{i,j};
    end
    if dz
        wfe = pupilOffset(wfe,Sensitivity.nomMask,dz/scale);
    end
    wfe(~Sensitivity.nomMask) = NaN;
    img=imagesc(scale*wfe); axis image;
    h = colorbar;
    title(h,Sensitivity.unitLabels{i});
%    set(img,'AlphaData',Sensitivity.nomMask);
    if dz
        title([Sensitivity.labels{i} ' (offset ' num2str(dz) ')']);
    else
        title([Sensitivity.labels{i}]);
    end
end

savepng([label ' consolidated segment sensitivites']);

%% how do you compute how beam launchers moved?
% hard to notice with translations -- too small for numerical precision.
% need micron perturbations to really see them.
% i = perturbation, j = segment
% i = 1;
% j = 1;
% xyz = Geometry.m1BeamLaunchers{j};
% moved = transform3D(xyz,
% moved = Sensitivity.transform{i,j} * [xyz'; ones(1,size(xyz,1))];
% moved = moved(1:3,:)';
% delta = moved-xyz;
% disp(delta);
