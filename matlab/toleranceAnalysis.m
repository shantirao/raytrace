function [Sensitivity,surfaces] = toleranceAnalysis(label,Geometry,nPoints,showPlots)

if nargin < 5
    showPlots = true;
end
if nargin < 4
    nPoints = 127;
else
    nPoints = max(ceil(nPoints/2)-1,1);
end

%%
aperture = Geometry.aperture;
surfaces = Geometry.prescription;

%% raytrace and showPlots
source = sourceColumn(aperture,nPoints,1);

options.negative=true;
trace = raytrace(source,surfaces,options);
trace = trace(2:end); %skip the source rays
nomRays = trace{end};
[nomOPL,nomMask] = pupilOPL(nomRays,false);
[rmswfe,piston] = rmsWFE(nomOPL,nomMask);
%%
if showPlots
    figure(3); clf;
    N = numel(trace);
    for i=1:N
        subplot(2,N,i);
        plotSpot(trace{i});
        plotApertures(trace{i}.surface);axis image;
        if isfield(trace{i}.surface,'name')
            title(trace{i}.surface.name);
        end

        subplot(2,N,N+i);
        displayOPL(trace{i});
    end
    savepng([label ' spot diagrams']);
end
%% perturb each segment 1 microradian or 1 micrometer
% raytrace unit system is millimeters and radians, so translations and
% rotation have different perturbation amounts. sorry.
% the "true" in perturb(s,x,true) means use the local coordinates of the mirror segment
perturbationScale = 1e-3*[1,1,1,.001,.001,.001];
perturbations =diag(perturbationScale); % 1 nm and 1 nanoradians
labels = {'Local \Deltax','Local \Deltay','Local \Deltaz','Local \thetaX','Local \thetaY','Local \thetaZ'};
fileLabels = {'Local x 1um','Local y 1um','Local z 1um','Local thetaX 1urad','Local thetaY 1urad','Local roll 1urad'};
unitLabels = {'um/um','um/um','um/um','um/urad','um/urad','um/urad'};

fig0 = 7;
maxI = 6; % number of degrees of freedom
maxJ = numel(surfaces{1}.segments);

plotM = ceil(sqrt(maxJ));
plotN = ceil(maxJ/plotM);

[scale, units] = displayScaleFactor(nomRays);

if showPlots
    figure(4);clf;
    x = scale*pupilOffset(nomOPL,nomMask,-piston);
    x(~nomMask(:)) = NaN; % AlphaData only works on Matlab
    img = imagesc(x);axis image;
%    set(img,'AlphaData',nomMask);
    title(sprintf('%s%5.5g %s rms','Nominal to pupil relay ' ,scale*rmswfe,units));
    title(colorbar,units)
    savepng([label ' nominal WFE']);
end
%% Boilerplate and setup
Sensitivity = struct();
Sensitivity.perturbations = perturbations;
Sensitivity.labels = labels;
Sensitivity.prescription = surfaces;
Sensitivity.source = source;
Sensitivity.opl = nomOPL;
Sensitivity.N = maxJ;
Sensitivity.nomMask = nomMask;
Sensitivity.nomrmsWFE = std(nomOPL(nomMask(:)));
Sensitivity.labels = labels;
Sensitivity.unitLabels = unitLabels;
Sensitivity.scale = scale;
Sensitivity.displayUnits = units;
%% M1 segments. Replace one at a time with a perturbed version
mask = cell(maxI,maxJ);
wfe = cell(maxI,maxJ);
rmswfe = cell(maxI,maxJ);
transform = cell(maxI,maxJ);

disp('this could take a while ...')
tic
% this is not a particularly fast way to do parallel processing because Octave takes a long time to start.
if exist('OCTAVE_VERSION', 'builtin') ~= 0
  if ~exist('pl',"dir")
    mkdir('pl');
  endif
  cmds = {};
  save("-z","pl/0","source");
  for i=1:maxI
    for j=1:maxJ
      sp = surfaces;
      delta = perturbations(i,:);
      s = surfaces{1}.segments{j};
      s = perturb(s,delta);
      sp{1}.segments{j} = s;
      pdata = sprintf('pl/%d_%d',i,j);
      save("-binary",pdata,"sp","options");
%      cmds{end+1} = ['octave --quiet --eval "load(\"pl0\");for j=1:' num2str(maxJ) ', load([\"pl/in_' num2str(i) '_\" num2str(j)]); prays=lastCell(raytrace(source,sp,options)); [opl,m] = pupilOPL(prays,false); save(\"-z\",[\"pl/out_' num2str(i) '_\" num2str(j)],\"opl\",\"m\");end"'];
      cmds{end+1} = ['load("pl/0"); load("' pdata '"); prays=lastCell(raytrace(source,sp,options)); [opl,m] = pupilOPL(prays,false); save("-z","' pdata '","opl","m");'];
    end
  end
  parallelRun(cmds);
  for i=1:maxI
    for j=1:maxJ
      pdata = sprintf('pl/%d_%d',i,j);
      load(pdata); % opl, m
      m = and(nomMask,m); %boolean and
      dOPL = pupilOffset(opl,m,nomOPL,@minus,0);
      rmswfe{i,j} = std(dOPL(m(:)));
      mask{i,j} = m;
      wfe{i,j} = dOPL;
      transform{i,j} = s.transform;
    end
  end
else
  for i=1:maxI
      parfor j=1:maxJ
          sp = surfaces;
          delta = perturbations(i,:);
          s = surfaces{1}.segments{j};
          s = perturb(s,delta);
          sp{1}.segments{j} = s;
          prays=lastCell(raytrace(source,sp,options));
          [opl,m] = pupilOPL(prays,false);
          m = and(nomMask,m); %boolean and
          dOPL = pupilOffset(opl,m,nomOPL,@minus,0);
          rmswfe{i,j} = std(dOPL(m(:)));
          mask{i,j} = m;
          wfe{i,j} = dOPL;
          transform{i,j} = s.transform;
      end
  end
end
toc
disp('... done calculating sensitivities')
Sensitivity.wfe = wfe;
Sensitivity.mask = mask;
Sensitivity.rmswfe = rmswfe;
Sensitivity.transform = transform;

%% M1 showPlots
if showPlots
    for i=1:maxI
        figure(fig0+i-1);clf;
        for j=1:maxJ
            subplot(plotM,plotN,j); %maxJ,maxI,(j-1)*maxI+i
            x = scale*wfe{i,j};
            x(~mask{i,j}(:)) = NaN;
            img = imagesc(x);axis image;
%            c = colorbar;
%            set(img,'AlphaData',mask{i,j});
            title(sprintf('%s %.4g %s rms/nm',surfaces{1}.segments{j}.name,scale*rmswfe{i,j},units));
        end
        if exist('majortitle'), majortitle([labels{i} ' [' unitLabels{i} ']'], 'fontsize',14);end
        savepng([label ' segment dWFE ' fileLabels{i}]);
    end
    if exist('tilefigs'), tilefigs(fig0:fig0+maxI-1); end
end

%% Consolidated
if showPlots
    figure(fig0+maxI+1);clf;
    offsets = [0 0 2 0 0 0];
    % something is wrong with j=14, i=3

    for i=1:6
        subplot(2,3,i);
        wfe = Sensitivity.wfe{i,1};
        dz = offsets(i);
        for j=2:maxJ
            wfe = wfe + Sensitivity.wfe{i,j};
        end
        if dz
            wfe = pupilOffset(wfe,nomMask,dz/scale);
        end
        x = scale*wfe;
        x(nomMask(:))=NaN;
        img=imagesc(x); axis image;
        h = colorbar;
        title(h,unitLabels{i});
%        set(img,'AlphaData',nomMask);
        if dz
            title([labels{i} ' (offset ' num2str(dz) ')']);
        else
            title([labels{i}]);
        end
    end

    savepng([label ' consolidated segment sensitivites']);
end

%% M2
mask = cell(maxI,1);
wfe = cell(maxI,1);
rmswfe = cell(maxI,1);
transform = cell(maxI,1);

tic
for i=1:maxI
    sp = surfaces;
    s = perturb(surfaces{2},perturbations(i,:),true);
    sp{2} = s;
    [prays]=lastCell(raytrace(source,sp,options));
    [opl,m] = pupilOPL(prays,false);
    m = and(nomMask,m); %boolean and
    dOPL = pupilOffset(opl,m,nomOPL,@minus,0);
    mask{i} = m;
    wfe{i} = dOPL;
    rmswfe{i} = std(dOPL(m(:)));
    transform{i} = s.transform;
end
toc

Sensitivity.m2wfe = wfe;
Sensitivity.m2mask = mask;
Sensitivity.m2rmswfe = rmswfe;
Sensitivity.m2transform = transform;

%% M2 showPlots
if showPlots
    figure(fig0+maxI+2);clf;
    for i=1:maxI
        subplot(2,3,i); %maxJ,maxI,(j-1)*maxI+i
        x = scale*wfe{i};
        x(mask{i}(:)) = NaN;
        img = imagesc(x);axis image;
%        c = colorbar;
%        set(img,'AlphaData',mask{i});
        title(sprintf('%s %.4g %s rms/nm',labels{i},scale*rmswfe{i},units));
    end
    if exist('majortitle'), majortitle(surfaces{2}.name , 'fontsize',14); end
    savepng([label ' ' surfaces{2}.name ' dWFE ']);
end

end
