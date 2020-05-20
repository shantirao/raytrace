function Q = rotationMatrixCq(axis,cq)
%rotationMatrixC a 3x3 Euler rotation matrix that operates on column vectors
% assuming you already have the cosine of the angle
% B = (rotationMatrix(axis,cos(angle)) * A')'
%
% Input:   axis [1x3]: Rotation axis (normalized)
%          angle  [1]: Rotation angle (radians)
%
% Output:   Q [3x3]: Euler rotation matrix
% -----------------------------------------------------------------------------
if abs(cq)>1e-18,
%     cq   = cos(angle); % cos(theta)
    sq   = sqrt(cq*cq); % sin(theta)
    omcq = 1-cq;      % 1 - cos(theta)

    Q = [ omcq*axis(1)^2+cq, ...
        omcq*axis(1)*axis(2)-sq*axis(3), ...
        omcq*axis(1)*axis(3)+sq*axis(2); ...
        omcq*axis(2)*axis(1)+sq*axis(3), ...
        omcq*axis(2)^2+cq, ...
        omcq*axis(2)*axis(3)-sq*axis(1); ...
        omcq*axis(3)*axis(1)-sq*axis(2), ...
        omcq*axis(3)*axis(2)+sq*axis(1), ...
        omcq*axis(3)^2+cq];
else
    Q = eye(3);
end

end 

