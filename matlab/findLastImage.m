function index = findLastImage(surfaces,sourceAperture,display)
%% find an exit pupil that is an optical conjugate of the first mirror
% inputs: surfaces is a cell array of optical surfaces, sourceAperture is a
% structure describing the ray origins
% outputs: pupil is the surface, index is where it goes in the surface 
% list, h is a plot handle. 

if nargin < 3
    display = false;
end

%% setup 
src = columnSource(sourceAperture,0,1); % single ray at center of aperture
options.segments = false;
options.aperture = false;
%% find an image by tracing from the source aperture. 

scale = src.aperture/100;
ri = src;
ri.N = 2;
ri.direction = repmat(src.direction,2,1);
ri.position = [src.position;src.position+max(scale)*src.local(1,:)]; % 1mm offset
ri.valid = [true;true];
ri.opl = [0;0];

[r,~,~,~,t1] = raytrace(ri,surfaces,options);

if display
    if isnumeric(display)
        figure(display); clf;
    elseif ~gcf
        figure();
    end
    
    subplot(2,5,[1,6]); 
    plotRays(t1,'b');
    plotApertures(surfaces,true);
    axis equal;
end

% find the last image in prescription -- surface where points are closest
% together
err = 1;
index = numel(surfaces)+1;
while index > 2 && err > 1e-9
% are they parallel?
    r1 = t1{index};
    index = index - 1;
    err = sum(diff(r1.position,1).^2);
    %abs(1-dot(r1.direction(1,:),r1.direction(2,:)));
end

% image position
[pimg,d] = lineIntersection(r1.position,r1.direction);
disp(sprintf('Image after surface %d(%s) at %d %d %d',index,r1.surface.name,pimg));
if display
    hold on; 
    scatter3(pimg(1),pimg(2),pimg(3),'r');  
    text(pimg(1),pimg(2),pimg(3),'Focus point','color','r','HorizontalAlignment','center');  
    hold off;
end