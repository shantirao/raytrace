%% Design parameters. Dimensions in mm.
apertureDiameter = 6500; 
Fnum = 11;
secondaryDistance = 4000;
backFocalLength = 4000; % distance behind m1 to prime focus, [mm]
segmentSize = 1500;
segmentSpacing = 40;
% mag = 29.17, not used
[prescription,aperture] = makeRCTelescope(apertureDiameter,Fnum,secondaryDistance,backFocalLength);
prescription{1}.segments = makeHexSegments(prescription{1},[0,0],segmentSize,segmentSpacing,pi/6,3);
aperture.units = 'mm';
label = 'HabeEX RC';
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
save Geometry;
save surfaceCenters;
%%
disp('M1 beam launchers')
for i=1:numel(Geometry.prescription{1}.segments)
    disp(Geometry.prescription{1}.segments{i}.name)
    disp(sprintf('[%g, %g, %g]\n',Geometry.m1BeamLaunchers{i}'))
end
disp('M2 Corner cubes')
    disp(sprintf('[%g, %g, %g]\n',Geometry.m2CornerCubes'))
    
%%
nPoints = 1024;
tic
[Sensitivity,surfaces] = toleranceAnalysis(label,Geometry,nPoints,true);
toc
%%
save([label 'Sensitivity'],'Sensitivity');

%% 
figure(2);clf;
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
    img=imagesc(scale*wfe); axis image;
    h = colorbar;
    title(h,Sensitivity.unitLabels{i});
    set(img,'AlphaData',Sensitivity.nomMask);
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