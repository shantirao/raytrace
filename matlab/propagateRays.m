function [rays, distance, projection, height] = propagateRays(rays, surface, options)
% see raytrace.m for instructions

optAperture = true;
optNegRays = false;
optSegments = true;
optParallel = false;
optPrecision = 1e-12; %< 1 fm. numerical: builtin('eps') * 2;
optDebug = false;
optFlatDeformations = false;

if nargin > 2
    if isfield(options,'aperture'), optAperture = options.aperture; end
    if isfield(options,'negative'), optNegRays = options.negative; end
    if isfield(options,'segments'), optSegments = options.segments; end
    if isfield(options,'parallel'), optParallel = options.parallel; end
    if isfield(options,'precision'), optPrecision = options.precision; end
    if isfield(options,'debug'), optDebug = options.debug; end
    if isfield(options,'flatDeformation'), optFlatDeformations = options.flatDeformation; end
end

maxIterations = 12;
nRays = size(rays.position,1);
% note that rays.direction is scaled by index^2
direction = (rays.direction); % could be a gpuArray
if ~isfield(rays,'status'), rays.status = ones(nRays,1); end
if ~isfield(rays,'valid'), rays.valid = true(nRays,1); end
if isfield(rays,'chief'), chief = rays.chief; else chief = 1; end;
if isfield(rays,'n2'), n2 = rays.n2; else n2=1; end
if isfield(rays,'polarization')
    polarization = [ones(1,nRays) zeros(nRays,3)] ;
end

label = '';
if isfield(surface,'label'), label = surface.label; end
if isfield(surface,'name'), label = surface.name; end
if ~isfield(surface,'type'), surface.type = 'stop'; end
if strcmp(surface.type ,'virtual')
  optNegRays = true;
end

if optDebug
  disp(sprintf('%s %s @ %d %d %d -> %d %d %d', label, surface.type, surface.position, surface.direction));
  disp(sprintf('  o  %.9d %.9d %.9d -> %.9d %.9d %.9d',rays.position(chief,:),rays.direction(chief,:)))
end

status = rays.status;
valid =  rays.valid; % could be a gpuArray
valid(chief) = true; % always trace the chief ray, and reset its valid status later. Only gets dropped from central obscurations

% calculate ray lengths L, error, new positions, surface normals N
% isDistorted = isfield(surface,'zernike') || isfield(surface,'grid') || isfield(surface,'mesh');
if isfield(surface,'curvature') && ~isfield(surface,'cuy')
    surface.cuy = surface.curvature;
elseif isfield(surface,'radius') && ~isfield(surface,'cuy') && surface.radius ~= 0
    surface.cuy = 1/surface.radius;
else
    surface.cuy = 0;
end
isFlat = ~isfield(surface,'cuy') || abs(surface.cuy) < optPrecision;
isSphere = ~isFlat && (~isfield(surface,'K') || abs(surface.K) < optPrecision);
isASphere = isfield(surface,'asphere') ;
isDistorted = isfield(surface,'deformation');
isZernike = isfield(surface,'zernike');
isNotSimple = isASphere || isDistorted || isZernike;

% in this branch, N is the surface normal vector. For planes it's a
% single vector. For curved surfaces, it's Nx3, a 3-vector for each
% intersection point.
% Bill Breckinridge at JPL used the equations to find the exact intersectio
% of a ray with a conic surface for the MACOS program, which was too much
% trouble to understand in FORTRAN, so I started over.
% Aspheres and other distortions are handled by iteratively transferring
% the surface intersection to a local tangent plane. It's arduous to do it
% in closed form, so you don't often see it derived in optics textbooks.
% It's just geometry. Iteration proceeds until the precision target
% (default 1e-14 of whatever units you're using) is reached.

%% find the Ray-Surface Intersection (RSI)
% note: N is the local surface normal at the point of intersection, with
% length 1. The ray directions are the length of the index of refraction.
% The ray bundle also carries around "n2", or "index squared" for the
% current medium.
% surface.index is the index of refraction of the next medium
position = rays.position - surface.position;

if isFlat && ~isNotSimple % simple plane
   N = surface.direction;
   bo = (N * direction')' ;
%   co = (N * bsxfun(@minus,surface.position,rays.position)')';
   co = -(N * (position)')';
   RSI = (abs(bo) > optPrecision);
   rays.status(valid & ~RSI) = -1; % No RSI
   valid = valid & RSI;
   L = co;
   L(RSI) = (L(RSI)./bo(RSI));  %avoid divide by zero
   P = projection = position + (direction.*L);
   height = 0;
%   newPosition = P + surface.position;
%        nearestNeighbor
   % dot products are a*b*cos(theta)
   % ray direction vectors are the length of the index of refraction
   % so NA is just the length of the component of the ray direction vector
   % orthogonal to the surface direction vector!
   % NA = N*sin(theta)
   nsint = sqrt(sum(multicross(direction,surface.direction).^2,2));
   if ~isempty(nsint)
    rays.NA = max(nsint(valid));
   end

else %not simple: distorted flat, or a conic, or an asphere
    % Initialize sign conventions
    P = position;
    Nv = surface.direction;

    if isFlat %
      c = 0;
      K = 0;
    else
      % The OSML framework wants the surface direction vector to point toward the active side of the surface
      % and for curvature to be positive when the direction points toward the center of the semi-latus-rectum.
      % However, CODE V and Zemax don't care which way the normal points. So if the surface direction and the rays
      % point in the same direction, we need to flip the logic around, and

      if isSphere
          K = 0;
      else
          K = surface.K;
      end

      c  = surface.cuy;
      sameSide = sign(P(chief,:) * Nv') == 1;
      sameDirection = sign(direction(chief,:) * Nv') == 1;
      if optDebug
        disp(sprintf('    side %d ddir %d',sameSide, sameDirection));
      end
      if sameDirection
        Nv = -Nv;
        c = -c;
      end
      if isfield(surface,'flip') && surface.flip % Nv is supposed to indicate the active surface
        if optDebug, disp('    flip'), end;
        Nv = -Nv;
        c = -c;
      end
    end

    if isSphere % K = 0 simplifies the equations
        Ao = -c*n2;
        Bo = (direction*Nv') - c*sum(P.*direction,2);
        Co = c*sum(P.*P,2) - 2*(P * Nv');
    else % K = anything; reduces to the previous case for K = 0 or a flat for cuy = 0
        RN   = direction*Nv';  % scalar ray . direction
        PN = (P*Nv');
        cKPN = c*K*PN;          % scalar position . opticdirection

        Ao = -c*(n2 + K*RN.^2);
        Bo = RN .* (1 - cKPN) - c*sum(P.*direction,2);
        Co = c*sum(P.*P,2) + PN.*(cKPN-2); % cK (P*Nv').^2 - 2(P*nV')
    end

    Wo = Bo.^2+Ao.*Co;

    % check if RSI exist
    RSI = (Wo >= 0);
    rays.status(valid & ~RSI) = -1; % No RSI
    valid = valid & RSI;

    Wo(RSI) = sqrt(Wo(RSI)); % avoid negative number with sqrt
    sw = sign(Bo);
    sw(sw == 0)=1;
    Qo = Bo + sw.*Wo;

    decision = Bo>=0;

    %suppress divide-by-zero
    if isSphere %Ao is a scalar
        t1 = RSI & decision;
        t1(t1) = abs(Ao) > optPrecision;

        t2 = RSI & ~decision;
        Qot2 = Qo(t2);
        t2(t2) = abs(Qot2) > optPrecision;

        L=zeros(nRays,1); % gpuArray.
        L(t1) = -Qo(t1)./Ao;
        L(t2) = Co(t2)./Qot2;
    else %Ao varies across the field
        t1 = RSI & decision;
        Aot1 = Ao(t1);
        t1(t1) = abs(Aot1) > optPrecision;

        t2 = RSI & ~decision;
        Qot2 = Qo(t2);
        t2(t2) = abs(Qot2) > optPrecision;

        L=zeros(nRays,1); %gpuArray.
        L(t1) = -Qo(t1)./Aot1;
        L(t2) = Co(t2)./Qot2;
    end
    if optDebug
      disp(sprintf('  L  %.9d',L(chief)));
    endif

%    disp(sprintf('%s %d %d %d %d %d',surface.name,L));
    P = position + (direction.*L);

% c is usually small, so multiplying it by P early removes precision in the
% ray position. P is multiplied by a unitary vector (Nv), which also can
% have 1 or 2 elements that are close to zero. Multiplying Nv*P first,
% and then scaling by c, preserves more precision.
    if ~isNotSimple
      if isFlat  % [1x3], same for every ray
          N = Nv;
      elseif isSphere % [M x 3]
          N = Nv - c*P; %bsxfun(@minus,Nv,c*P);
      else % [M x 3]
          N = Nv.*((1-(Nv*P')'*K*c))-c*P; %bsxfun(@times,Nv,(1-c*K*(Nv*P')'))-c*P;
      end
      %not an asphere. Easy to calculate these
      N = N./sqrt(sum(N.^2,2));
      height = (Nv * P')';
      projection = P - height .* Nv; %bsxfun(@times,height,surface.direction); % in global coordinates
   else
    %asphere terms or deformation; iterative solution for both height
    %and the surface angle. You could do asphere=[x] to see if it matches
    %a parabola exactly.
    %do this again because now we use the slower sag equation, and don't flip the vertex normal around
        c = surface.cuy;
        Nv = surface.direction;
%        oldP = position - surface.position;
        c2 = c*c;
        k1c2=(K + 1)*c2;
        Pa = P(RSI,:);
        pos = position(RSI,:);
        ddir = direction(RSI,:);
        if c == 0 % plane
           z = zeros(size(Pa,1),1); % 0 * height;
        else % sqrt requires numbers that actually work -- nonvalid elements have weird values that force conversion to complex numbers
           z = (Nv * Pa')';
        end
        for j=1:maxIterations
%          Prev = Pa;
            height = (Nv * Pa')';
            projection = Pa - height.*Nv; %

            if isfield(surface,'local')
                xy = projection * surface.local(1:2,:)';
            else
                error('a surface needs a local coordinate system to evaluate Zernike coefficients');
            end
            if isfield(surface,'center') % technically, a symmetric asphere doesn't need uv
                uv = xy - surface.center(1:2);
            else
                uv = xy;
            end

%            if isfield(surface,'center') % technically, a symmetric asphere doesn't need uv
%                uv = projection * surface.local(:,1:2)' - surface.center(1:2);
%            elseif isfield(surface,'local')
%                uv = projection * surface.local(:,1:2)';
%            else
%                error('a surface needs a local coordinate system to evaluate Zernike coefficients');
%            end

            s2 = sum(xy.^2,2); % dot(projection,projection)
            if ~isFlat
                z = c*s2./(1+sqrt(1-(k1c2*s2))); %sag equation
                dudv = c*xy./sqrt(1-k1c2*s2); %directional gradient
            else
                z = zeros(size(Pa,1),1);
                dudv = zeros(size(xy)); % no slope since
            end

            % iterative corrections: find the intersection of the ray with a plane on the calculated surface

            if isASphere
                %aspheric component of curvature, for adding to the normal
                %vector  z(s) = a2 s^2 + a4 s^4 + a6 s^6 + ...
                %    dz/ds(s) = s * (2*a2 + 4*a4*s^2 + 6*a6*s^4)
                %         (one factor of s missing for now. Save it for
                %         later, and it provides direction
                f = 2;
                radialDerivative = zeros(nRays,1);
                for a = surface.asphere
                  radialDerivative = radialDerivative + f * a * sn;
                  sn = sn .* s2;
                  dz = dz + a * sn;
                  f = f + 2;
                end
                Pnear = projection + z * Nv;% P is the nearest point on the parent conic

                z += dz;
                dxdy = radialDerivative .* xy ./ sqrt(s2); %x and y components
                dudv += dxdy;

                % Nearest point on the asphere to the last estimate
%                Pnear = tangentProjection + bsxfun(@times,(dz+z),Nv);
%                N = N - bsxfun(@times,(Nv * N')' .* radialDerivative,tangentProjection);
               % update the normal vector.
                % add the radial aspheric component of the derivative to N
                % scale the aspheric component to match the conic normal magnitude
                % This may appear to be bonkers. For a surface defined by the locus of points P(x,y) = [x, y, z(x,y)], the surface
                % normal direction is always [dz/dx, dz/dy, 1]. You can derive it from N = cross(dP/dx, dP/dy). With linear
                % perturbations, P(x,y) = P_conic(x,y) + P_asphere(x,y) we get N(x,y) = N_conic(x,y) + N_asphere(x,y). That's why
                % it's important not to normalize N until all the perturbations are calculated.
            end

            % weird deformations. If "center" is specified, deformations are relative to the aperture, not the parent surface
            % if not, define a local coordinate system
            if isZernike || isDistorted

            if isZernike % sym
              if isfield(surface.zernike,'center')
                uv = xy - surface.zernike.center;
              end

              if isfield(surface.zernike,'Noll')
                [dz dxdy] = zernikeIndex(surface.zernike.Noll, uv ,surface.zernike.radius);
%                  ap = surface.zernike.radius;
%                  if isfield(surface,'zernikeRadius')
%                    ap = surface.zernikeRadius;
%                  elseif isnumeric(surface.aperture)
%                      ap = max(surface.aperture); %circle or annulus
%                  else
%                      ap = surface.aperture.radius(end); %new form- structure describes aperture. much better
%                  end
%                  error('Noll zernikes not implemented. See zernikeIndex.m ')
              elseif strcmpi(surface.zernike.type,'XY')
                [dz dxdy] = zernikeXY(surface.zernike, uv, surface.zernike.R);
              else
                error('unknown type of zernikes deformation')
              end

                % surface.local: [2x3]
                % N: [Nx3]
                % dxdy: [Nx2}
              z += dz;
              dudv += dxdy;
%                N = N - (dxdy*surface.local(1:2,:)); % normal is perpendicular to directional slope.
% z = f(x,y)  -> The gradient of z - f(x,y) = 0 is the normal vector. f(x,y) = sag(x,y) + zernike(x,y), so add the derivatives from sag (dudv) and zernike (dxdy)
            end

            if isDistorted % Octave has different functions for triangulated meshes https://docs.octave.org/v10.1.0/Identifying-Points-in-Triangulation.html
                % interpolate a triangulated mesh, where Pnear moves
                % with the interpolated average of the distortion
                % points, and N moves with the interpolated difference

                % vertex indices and barycenters, in local coordinates

                % have to iterate over the mesh because a point can be
                % in more than one triangle if it's at a vertex
                for i=1:size(N,1)
                    [triangle,barycenter] = pointLocation(surface.deformation.triangulation,uv(i,:));
                    triVals = surface.deformation.displacement(surface.deformation.triangulation(triangle,:),:);
                    dz(i) += barycenter * triVals; %(bc',triVals')';
%                    dP = barycenter * triVals; %(bc',triVals')';
                    % dP = mean(dP,1); %multiple triangles get averaged
                    % translates the surface intersection point.
                    % pointLocation used to return more than one triangle

                    % if deformations are surface-normal, instead of
                    % tangent-normal, you would do dP * normr(N) here
%                    Pnear(i,:) = Pnear(i,:) + dP * Nv;

                    % old: the surface normal moves a bit, too. It can be
                    % expressed as a rotation. Use the curl of the
                    % triangular facets to rotate the normal vector
%                     rotvec = mean(surface.deformation.curl(vi,:),1);
%                     Q = rotationMatrix(normr(rotvec),norm(rotvec));
%                     N(i,:) = Q * N(i,:)';

                if ~optFlatDeformations % change WFE without steering the beams
                    % use flatDeformations = true if tangent-normal
                    % deformations are small, and not gridded finely enough
                    % to calculate the beamwalk

                    % if a point falls on an edge or vertex, the path length is
                    % correct, but the tilt comes from only one of the triangles associated
                    % with that vertex. If any of the barycenters are less than a precision
                    % threshold, average the tilts from the adjacent triangles as well.
                    notIt = barycenter < optPrecision;
                    s = sum(notIt);
                    if s == 2 %vertex - single
                        vertex = surface.deformation.triangulation.ConnectivityList(triangle,~notIt);
                        dxdy(i,:) += surface.deformation.vertexSlope(vertex,:);
                    elseif s == 1 %edge - 2 vertices are close
                        vertex = surface.deformation.triangulation.ConnectivityList(triangle,~notIt);
                        edges = all(surface.deformation.edges == sort(vertex),2);
                        %weighted average of vertices? No, better to use
                        %the average of the adjacent triangles
                        %dxdy = sum(surface.deformation.edgeSlope(vertex,:).*barycenter(~notIt),1);
                        dxdy(i,:) += surface.deformation.edgeSlope(edges,:);
                    else % look up the slope for that particular vertex
                        dxdy(i,:) += surface.deformation.slope(triangle,:);
                    end
%                    N(i,:) = N(i,:) - dxdy * surface.local;
                end
                end
            end
            end
            Pnear = projection + z .* Nv ;
            Nnear = Nv - (dudv)*surface.local(1:2,:); %  (-df/dx, -df/dy, +1) in global coordinates
            %find the ray-surface intersection for Pnear and Nnear
            Lnew = sqrt(sum((Pnear - pos).^2,2));
%            errest = abs(P(valid,:) - Pnear(valid,:));
%project the ray onto the local tangent plane. Is it the point we found?
            La = sum((Pnear - pos) .* Nnear,2) ./ sum(ddir .* Nnear,2); % don't have to normalize Nnear because it's in the denominator
            Pa = pos + ddir .* La ; % starting point for next iteration -- has to be a point on the original ray
%errest = max(abs(pos + ddir.*Lnew - Pnear));
%            errest = abs(Prev(:) - Pa(:));
            errest = abs(Pnear(:) - Pa(:));
            if all(errest(:) < optPrecision) % errest should be identically 0 for on-axis rays
              break
            end
        end
        if nargout > 3
              height = (Nv * Pa')';
        end
        P(RSI,:) = Pa;
        L(RSI,:) = La;
        N = zeros(size(position));
        N(RSI,:) = Nnear./sqrt(sum(Nnear.^2,2)); %normalize
    end

%    delta = newPosition - position;
%    L = sqrt(sum(delta.^2,2));
%    N = N./sqrt(sum(N.^2,2)); %normalize

if optDebug
  if isZernike && strcmp(surface.zernike.type,'XY')
%    verify the location matches the sag equation here
    disp(surface.name)
    height = (Nv * P')';
    projection = P - height .* Nv;
    uv2 = projection * surface.local';
    s3 = sum(uv2.^2,2);
    h2 = c*s3./(1+sqrt(1-(k1c2*s3)));
    dz = zernikeXY(surface.zernike, uv, surface.zernike.R);
    z2 = h2 + dz;
    P2 = uv2 * surface.local + (z2) * Nv;
    disp(P2 - P)
  end
end


%    N(sameDirection(:),:) = - N(sameDirection(:),:);
    % possible here because the rays might point the wrong way
end

% projection is in global coordinates
% center is in tangent plane coordinates
% uv = projection * surface.local' gives local coordinates
% uv = (projection*surface.local') - surface.center gives aperture coordinates
%%
%if isfield(surface,'cuy') && surface.cuy < 0
%    height = -height;
%end

% inAP = apertureTest(surface, projection)
%% Aperture test, in tangent plane space. Not optimized for
% computational efficiency because there are some extra
% translations into and out of the local coordinate frame
% center is in the uv plane -- if there's a third component, it's height
% any aperture other than a circle or an annulus needs a local coordinate system
inAP = 1;
uv = [];
if optAperture && (isfield(surface,'aperture') && ~isempty(surface.aperture) ||isfield(surface,'diameter'))
 [inAP, uv] = apertureTest(surface, projection);
elseif isfield(surface,'pixelSize')
    if isfield(surface,'local')
      if isfield(surface,'center')
        uv = projection * surface.local(1:2,:)' - surface.center(1:2);
      else
        uv = projection * surface.local(1:2,:)';
      end
    else
        uv = projection;
    end
end %uv only needed for pixel size

%% Map onto pixels
if isfield(surface,'pixelSize')
    rays.pixels(valid,:) = round(fix(2*(uv(valid,:))./surface.pixelSize)/2);
    rays.pixelSize = surface.pixelSize;
end

%% Found the RSI! Now what? Update the path length
%position = newPosition ;
%distance = n2 .* L;
%     opl = opl + distance;

%% Reflect or refract?
surfType = 0; % stop or reference
if isfield(surface,'type')
    surfType = surface.type;
    if strcmp(surfType,'reflect')
        surfType = 1;
    elseif strcmp(surfType,'refract')
        surfType = 2;
    end
end

if surfType == 2 %maybe it changes?
    if isfield(surface,'material')
        if ~isfield(rays,'wavelength')
            error('entering a dispersive medium. rays.wavelength must be defined');
        end
        if ~isfield(rays,'units')
            error('please define units for the wavelength, such as rays.units = ''mm''');
        end
        ne2 = surface.material.evaluate(surface.material, rays.wavelength, rays.units);
    elseif isfield(surface,'index')
        ne2 = surface.index*surface.index;
    else
        error('refractive surfaces must have an index of refraction');
    end
    rays.n2 = ne2;
else
    ne2 = n2; % no change, but maybe the reflection grating needs it
end

if isfield(surface,'grating')
    if ~isfield(rays,'wavelength')
        error('must set rays.wavelength for gratings to work');
    end
    % for a grating defined by a planar surface, surface.grating is a
    % vector of length (order / spacing), pointing in
    % the direction of the parallel planes.
    % If you can define the grating vector r3 as a function of position on
    % the tangent plane, then you can define an arbitrary hologram.
    % See C J Mitchell 1981 J. Opt. 12 301
    % "Generalized ray-tracing for diffraction gratings of arbitrary form"
    % variable names changed
    % N          n0       [Nx3] surface-normal (length 1)
    % direction  mu1*n1   [Nx3] ray vector (length index-of-refraction)
    % grating    m*n3/g   [1x3] order * grating vector / grating period
    % r2         mu2*n2   [Nx3] output ray vector (length index of the next medium)
    % C          m*lambda/g     order * wavelength / grating period
    lambda = rays.wavelength;
    grating = surface.grating;
    if isfield(surface,'order'), grating = grating * surface.order; end

%     mu2 = 1;
%     if isfield(surface,'n'), mu2 = surface.n; end

    % note: dot(N,direction) is sum(N.*direction,2)
    n0r1 = sum(N.*direction,2);
    r3 = lambda.*grating;
    n0r3 = sum(N.*r3,2);
    r1r3 = sum(direction.*r3,2); % each ray could have a different wavelength
    r3r3 = sum(r3.*r3,2); % order*lambda/pitch

    sgn = sign(n0r1); %refract
%     if isfield(surface,'n')
%         ne2 = surface.n*surface.n;
%     elseif surfType == 1 %reflect
%         ne2 = n2;
%     else
%         error('refractive gratings must have surface.n = index of refraction');
%     end

    if surfType == 1, sgn=-sgn; end %reflect;

    TIR = (n0r1-n0r3).^2 + ne2 - n2 -r3r3 - 2*r1r3;
    direction = direction - r3 + N.*(sgn.*sqrt(TIR) - n0r1 - n0r3);

    %whew!


else %reflect or refract? Note N, the local surface normal, has length 1
    switch surfType
        case 1 %reflect
            if isFlat && ~isNotSimple
                direction = bsxfun(@minus,direction,2*(N*direction')'*N);
            else % dot product = sum(N.*direction,2)
                delta = 2*bsxfun(@times,N,sum(N.*direction,2));
                direction = bsxfun(@minus,direction,delta);
            end
        case 2 %refract
            np=sum(bsxfun(@times,N,direction),2);
            sgn = sign(np);
            sgn(sgn==0) = 1;
            TIR = np.^2 + ne2 - n2;

            %error detection
            errTIR = TIR<0;
%             lost = valid;
%             lost(valid) = errTIR;
            rays.status(valid & errTIR) = -1;
            valid = valid & ~errTIR;

%             t = bsxfun(@times,N,np-surface.*sqrt(TIR));
            t = N .* (np-sgn.*sqrt(TIR));
            direction = direction - t;
%            direction(valid,:) = direction(valid,:) - t(valid,:);

            % new surface
%             rays.n2 = ne2;

        case 4 %to do: polarization: reflect and possibly refract
            error('sorry, polarization surfaces are not yet implemented');
            % make a new polarizationReflection() function
            if isfield(rays,'local') %project the TE vector onto the reflection plane
                TE = rays.local;
            else
                TE = bsxfun(@cross,dirnorm,N);
            end
            n = sqrt(n2);
            dirnorm = direction / n;
            s = bsxfun(@cross,dirnorm,N); %ray direction before reflection [N x 3]
            s = normr(s); % problem at normal incidence
            cosphi = sum(s.*rays.local,2)/ n; % dot product, scaled down by the incident ray
            c2phi = 2 * cosphi.^2 - 1;
            %could be a sign error here -- have to check that with examples
            sinphi =sum(bsxfun(@cross,rays.local,s).*direction,2);
            s2phi = 2*cosphi.*sinphi;
    %             rotationMatrix = [1 0 0 0; 0 c2phi s2phi 0; 0 -s2phi c2phi 0; 0 0 0 1];
            if isfield(rays,'polarization')
                if isfield(surface,'NK')
                    n1cosincidence = bsxfun(@dot,newDirection,N)/n;
                    [rp, rs, m11, m12, m33, m34] =  fresnel(n1cosincidence,n,surface.NK);
                else % perfect reflector
                    m11 = 1; m12 = 0; m33 = -1; m34 = 0;
                end
            end
            stokes = rays.polarization;%for ease of debugging
            stokes = muellerRotate(c2phi,s2phi,stokes);
            stokes = muellerMul(m11, m12, m33, m34,stokes);
            rays.polarization = stokes;
            rays.local = s; % new TE vector
            if isfield(rays,'local')
                xdotN = sum(bsxfun(@times,N,rays.local),2);
                rays.local = N .*xdotN + bsxfun(@times,rays.local,-1 + xdotN);
            end
            if isfield(rays,'phase')
                rays.phase = -rays.phase;
            end
            if isfield(rays,'polarization')
                %http://people.physics.tamu.edu/mgao/calculator/my.jqplot.fresnel.html
            end
            if isfield(rays,'local')
                xdotN = sum(bsxfun(@times,N,rays.local),2);
                rays.local = N .*xdotN + bsxfun(@times,rays.local,-1 + xdotN);
            end

        otherwise % do nothing -- no change to direction

    end
end

% Check for rays that shouldn't have traced
negRays = (L<0);
rays.status(valid & negRays) = -2; % Negative ray
if ~optNegRays
    valid = valid & ~negRays;
end
if optAperture
    rays.status(valid & ~inAP) = 0; % aperture
    valid = valid & inAP;
end

distance = n2 .* L;
if isfield(surface,'display'), rays.display = surface.display; end
if isfield(rays,'opl'), rays.opl = rays.opl + distance; else rays.opl = distance; end

%
%todo: transform the coordinate system that the rays drag along with them
% u = new direction
% v = original direction
% If u . v < 1, they're in opposite directions, so need a sign
if isfield(rays,'TE')
    TE = rays.TE';
    if size(TE,2) == 1
        TE = repmat(TE,1,nRays);
    end
    sOld = 1/sqrt(rays.n2); % need unitary rotatioon matrices, and ray direciton
    sNew = 1/sqrt(ne2);     % is scaled by the index of refraction
    sgn = 1;
    if surfType == 1
        sOld = -sOld;
    end
    for i=1:nRays
        if valid(i)
            % direction = r * (old direction)
            r = rotationBetween(sOld*rays.direction(i,:)',sNew*direction(i,:)');
            if surfType == 1  %reflect
                TE(:,i) = -(r*(TE(:,i)));
            else
                TE(:,i) = (r*(TE(:,i)));
            end
        end
    end
    rays.TE = TE';
end

v = valid;
v(chief) = true;

% show where the position would be for rays that missed the aperture, but don't propagate them any further.
%rays.position(v(:),:) = P(v(:),:) + surface.position;
%rays.position(~v(:),:) = NaN;
%rays.direction(v(:),:) = direction(v(:),:);
%rays.direction(~v(:),:) = NaN;
rays.position = P + surface.position;
rays.direction = direction;
rays.status = status;
if ~rays.valid(chief) % if this was knocked out earlier, mark it as invalid, even though we compute its path for other reasons.
  valid(chief) = false;
end
rays.valid = valid;
rays.normal = N;
rays.surface = surface;
if optDebug
  disp(sprintf('  i  %.9d %.9d %.9d -> %.9d %.9d %.9d',rays.position(chief,:),rays.direction(chief,:)))
  projection(~valid,:)=0;
  rays.projection = projection; %in the global coordinate system, relative to the vertex position
  rays.distance = distance;
  rays.height = height;
end


function R = rotationBetween(u, v)
% see https://math.stackexchange.com/questions/432057/how-can-i-calculate-a-4-times-4-rotation-matrix-to-match-a-4d-direction-vector/2161406#2161406
% from https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/2161631#2161631
% returns R such that v = R * u, and R * cross(u,v) = cross(u,v)
% u and v are both column vectors
    h = v+u;
    S = eye(3) - 2 * h * (h')/(h'*h);
    R = S - 2 * v * (v' * S)/(v'*v);
end


end

