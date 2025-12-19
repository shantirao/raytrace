function surface = perturb( surface, delta, inLocalCoordinates)
%perturb moves a surface location, direction, and local coordinate system.
%delta can be a structure describing an vector deltas to position,
% direction, or {angle, center, and axis}. The surface.transform output is
% suitable for moving the mesh point locations
% delta.angle       [1]: Rotation angle, in radians
% delta.axis      [1x3]: Rotation axis (normalized, defaults to global X)
% delta.center    [1x3]: Center of rotation (defaults to surface vertex)
% delta.position  [1x3]: Offset of vertex
% delta.direction [1x3]: Vector displacement to direction
% delta.scale       [1]: Scale the angle or displacement
% delta           [1x6]: [dx, dy, dz, thetaX, thetaY, thetaZ] about surface
%                        (in Zemax, X, Y, Z, U, V, W,)
%                        local coordinates, which default to global x,y,z
% The local coordinate system (surface.local) and reference axes
% (surface.reference) are also rotated because they are in the global
% coordinate system.

% Tranformation matrix: see https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html

if nargin < 3
    inLocalCoordinates = false;
end

if iscell(surface)
  for i=1:numel(surface)
    surface{i} = perturb(surface{i},delta, inLocalCoordinates);
  endfor
elseif iscell(delta)
    for i=1:numel(delta)
        surface=perturb(surface,delta{i},inLocalCoordinates);
    end
else
    %center of perturbation in global coordinates but the offset
    %is carried forward in the tangent plane for child segments to adjust
    %according to their particular location. If the elevation is known,
    %apertureCenter accounts for that and the transformation is in the parent
    %coordinate system, but about the center point of the aperture.
    %reference coordinates might be offset and tilted, whereas local coordinates travel with the parent vertex for an off-axis aperture
    if inLocalCoordinates && isfield(surface,'referenceCenter')
      center = surface.referenceCenter;
    else
      center = surface.position + apertureCenter(surface);
    end
    label = '';
    %passthrough to subpaertures from parent
    if isnumeric(delta) && numel(delta) == 16 && length(delta) == 4
        T = delta;
    else
        if isnumeric(delta)
            d = delta(1:3); %displacement
            if numel(delta) < 6
                r = [0,0,0];
            else
                r = delta(4:6); %rotation
            end
    %         surface.position = surface.position + d(1:3);
    %         move w.r.t. the surface reference axes
            if inLocalCoordinates
                if isfield(surface,'reference')
                    reference = surface.reference;
                elseif isfield(surface,'local')
                    reference = [surface.local; surface.direction];
                else
                    error('a surface needs reference axes to be perturbed this way')
                end
                Q = rotationMatrix(normr(reference(1,:)),r(1)) * ...
                    rotationMatrix(normr(reference(2,:)),r(2)) * ...
                    rotationMatrix(normr(reference(3,:)),r(3));
                d = (d * reference)';
            else
                Q = rotationMatrix([1 0 0],r(1)) * rotationMatrix([0 1 0],r(2)) * rotationMatrix([0 0 1],r(3));
%                 d = d(4:6);
            end
        else
            if inLocalCoordinates
                error('perturbations in local coordinates only work when the perturbation is a numeric vector');
            end

            if isfield(delta,'scale')
                scale = delta.scale;
            else
                scale = 1;
            end
            if isfield(delta,'label')
                label = [' ' delta.label];
            end
            %rotation first
            if isfield(delta,'direction') % vector offset
                %surface.direction is always normalized
                x = normr(surface.direction + scale*delta.direction);
                Q = rotationMatrixCq(normr(cross(surface.direction,x)),dot(surface.direction,x));
    %             surface.direction = x;
            elseif isfield(delta,'angle')
                Q = rotationMatrix(normr(delta.axis),scale*delta.angle);
                if isfield(delta,'center')
    %             	oldCenter = center;
                    center = delta.center; %override the optic center position)
                end
                % for verifying the transformation matrices below
    %             sp = center + (Q * (surface.position - center)')';
    %             sd = (Q * surface.direction')';
            else
                Q = eye(3);
            end

            %now translate to move the origin
            if isfield(delta,'position') % translation
                d = scale*delta.position;
    %             sp = sp + d;
            else
                d = [0 0 0];
            end
        end

        T = [ Q d(:); 0 0 0 1]; %simultaneous translate and rotate. note: last row is skew & scaling. don't mess.
        C2 = [eye(3) -center'; 0 0 0 1]; % translate
        C1 = [eye(3) center'; 0 0 0 1]; %translate back
        T = C1 * T * C2;
    end

    if isfield(surface,'transform') %surface.transform remembers the history of moving about
        surface.transform = surface.transform * T;
    else
        surface.transform = T;
    end

    % how to transform a global position vector X?
    % pos = (surface.transform * [X 1]')(1:3)
    % dir = (surface.transform * [X 0]')(1:3)
    % surface.center and segment.aperture are in the local coordinate
    % system so they do not transform.
    o = ones(size(surface.position,1),1);
    z = zeros(size(surface.direction,1),1);
    y = (T * [surface.position o]');
    surface.position = y(1:3,:)';
    y = (T * [surface.direction z]');
    surface.direction = y(1:3,:)';

    % local coordinates
    if isfield(surface,'local')
        z = zeros(size(surface.local,1),1);
        y = (T * ([surface.local'; z']));
        surface.local = y(1:3,:)';
    end

    if isfield(surface,'label')
        surface.label = [surface.label label];
    end

    if isfield(surface,'referenceCenter') % move the points
%      surface.referenceCenter = surface.transform * [0 0 0;0 0 0;0 0 0;surface.referenceCenter];
%      surface.referenceCenter = surface.transform * [surface.referenceCenter 1]';
      pts = surface.referenceCenter;
      a = T * [pts'; ones(1,size(pts,1))];
      surface.referenceCenter = a(1:end-1,:)';
    end
    if isfield(surface,'fiducials') % move the geometric datums
      pts = surface.fiducials;
      a = T * [pts'; ones(1,size(pts,1))];
      surface.fiducials = a(1:end-1,:)';
    end

    % reference axes - change direction but don't translate
    if isfield(surface,'reference')
        y = (T * ([surface.reference'; zeros(1,size(surface.reference,1))])); %0 0 0
        surface.reference = y(1:3,:)';
    end

    if isfield(surface,'segments')
        for i=1:numel(surface.segments)
            surface.segments{i} = perturb(surface.segments{i},T);
        end
    end
end
