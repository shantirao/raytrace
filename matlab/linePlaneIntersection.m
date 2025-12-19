function [position, distance] = linePlaneIntersection(linePos,lineDir, planePos, planeDir)
%tired of writing this out every time1
   bo = (planeDir * lineDir')' ;
   co = -(planeDir * (linePos - planePos)')';
   distance = co./bo;
   position= linePos + distance.*lineDir;

