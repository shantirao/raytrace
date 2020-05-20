example_303
%shift the lens over by 3 mm. Should change the wavefront a bit.
trace2 = raytrace(source,{perturb(s1,[3,0,0,0,0,0]),s2,s3});
std(trace{end}.opl - trace2{end}.opl)
%% shift direction vector
surface = struct('position',[0,0,0],'direction',normr([.5 .5 0]))
delta=struct('direction',[.1,-.1,0])
s = perturb(surface,delta)
b = s.transform * [surface.direction 0]'
x = norm(b(1:3)'-s.direction)
% should be within numerical error

%% rotation about a point
surface = struct('position',[.5 .5 0],'direction',[1 0 0])
delta = struct('center',[0 1 0],'axis',[0 0 1],'angle',pi/4) % center of rotation
s = perturb(surface,delta)
b = s.transform * [surface.position 1]';
disp(sprintf('transform * position = \n\t [%3.3f %3.3f %3.3f]',(b(1:3))));

%% debug
acosd(dot(s.direction,surface.direction))