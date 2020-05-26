set(gcf,'Units','normal');
% set(gca,'Position',[.05 .05 .9 .9]);
axis equal;
axis tight;
view([1,0,0]);
daspect([1,1,1]);
camroll(180);
cameratoolbar('SetCoordSys','y');
cameratoolbar('setmode', 'orbit') ;
grid on