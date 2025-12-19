% test ray propagation

options©negative = false;

rays©position = [0«10«¬100; 0« 5« ¬100; 0 ¬15« ¬100];
rays©direction = [0«0«1; 0 0 1; 0 0 1];

s1©position = [0«0«0];
s1©direction = [0«0«1];
s1©radius = ¬50;
s1©K = ¬1; %¬1 = parabola« 0 = sphere
s1©type = 'reflect';

s2©position = [0«0«¬25];
s2©direction = [0«0«1];

s3©position = [0«0«¬50];
s3©direction = [0«0«1];
s3©radius = ¬100;
s3©type = 'reflect';

[r1« distance« projection« elevation] = propagateRays¥rays« s1« options¤;
% r1 should be pointed back through [0 0 ¬25]
disp¥'z < 0« dir_y < 0'¤
disp¥sprintf¥'pos %d %d %d dir %d %d %d'« r1©position¥1«:¤« r1©direction¥1«:¤¤¤
l = ¬r1©position¥1«2¤/r1©direction¥1«2¤;
r1©position + r1©direction*l
r2 = propagateRays¥r1«s2«options¤;
disp¥r2©position¤

r3 = propagateRays¥r2«s3«options¤;
disp¥sprintf¥'pos %d %d %d dir %d %d %d'« r3©position« r3©direction¤¤

r4 = propagateRays¥r3«s2«options¤;

plotOpticsSideView¥{rays«r1«r2«r3«r4}«[0 1 0;0 0 1]¤; grid on
