function [rays, distance, projection, elevation] = propagateRays(rays, surface, options)
% see raytrace.m for instructions

optAperture = true;
optNegRays = false;
optSegments = true;
optParallel = false;
optPrecision = 1e-14;
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

maxIterations = 10;
nRays = size(rays.position,1);
position = (rays.position); % could be a gpuArray
% note that rays.direction is scaled by index^2 
direction = (rays.direction); % could be a gpuArray
if ~isfield(rays,'status'), rays.status = ones(nRays,1); end
if ~isfield(rays,'valid'), rays.valid = true(nRays,1); end
if isfield(rays,'n2'), n2 = rays.n2; else n2=1; end
if isfield(rays,'polarization')
    polarization = [ones(1,nRays) zeros(nRays,3)] ;
end

label = '';    
if isfield(surface,'label'), label = surface.label; end
if isfield(surface,'name'), label = surface.name; end
status = rays.status;
valid =  (rays.valid); % could be a gpuArray

oldValid = valid;

% calculate ray lengths L, error, new positions, surface normals N
% isDistorted = isfield(surface,'zernike') || isfield(surface,'grid') || isfield(surface,'mesh');
if isfield(surface,'curvature') && ~isfield(surface,'cuy')
    surface.cuy = surface.curvature;
elseif isfield(surface,'radius') && ~isfield(surface,'cuy') && surface.radius ~= 0
    surface.cuy = 1/surface.radius;
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
% surface.n is the index of refraction of the next medium
if isFlat && ~isNotSimple % simple plane
   N = surface.direction; 
   bo = (N * direction')' ;
   co = (N * bsxfun(@minus,surface.position,rays.position)')';
   RSI = (abs(bo) > optPrecision); 
   valid = valid & RSI;
   L = co;
   L(valid) = (L(valid)./bo(valid));
   newPosition = position + bsxfun(@times,direction,L);
   elevation = 0;
   projection = bsxfun(@minus,newPosition,surface.position);
   if isfield(surface,'local')
      projection = projection * surface.local';
   end

%        nearestNeighbor
   % dot products are a*b*cos(theta)
   % ray direction vectors are the length of the index of refraction
   % so NA is just the length of the component of the ray direction vector 
   % orthogonal to the surface direction vector!
   % NA = N*sin(theta)
   nsint = sum(bsxfun(@minus,direction,surface.direction).^2,2);
   nsint = nsint(valid);
   if ~isempty(nsint)
    rays.NA = (max(nsint)-min(nsint))/2;
   end 

else %not simple: distorted flat, or a conic, or an asphere  
    % Initialize sign conventions
    if isFlat
        c = 0;
        K = 0;
        Outside = false; 
        Nv = surface.direction;
        sgn = -1;
    else
        % Outside means find the solution for the outside of the conic surface
        % using the sign is a terrible signal. Just set a parameter for
        % that like "isConvex"
        Outside = surface.cuy < 0; 
        if isSphere
            K = 0;
        else
            K = surface.K;
        end
        if Outside
            Nv = -surface.direction;
            c  = -surface.cuy;
            sgn = 1;
        elseif isfield(surface,'convex') && surface.convex
            % surface.direction points the direction of the conic, and the
            % radius of curvature is positive.
            Outside = true;
            Nv = surface.direction;
            c  = surface.cuy;
            sgn = 1;
        else
            Nv = surface.direction;
            c  = surface.cuy;
            sgn = -1;
        end            
    end

%     P = bsxfun(@minus,position,surface.position);
    P = position - surface.position;

    if isSphere % K = 0 simplifies the equations
        Ao = -c*n2;
        Bo = (direction*Nv') - c*sum(P.*direction,2);
        Co = c*sum(P.*P,2) - 2*(P * Nv');
    else % K = anything; reduces to the previous case for K = 0 or a flat for cuy = 0
        RN   = direction*Nv';  % scalar ray . direction
        %cP   = c*P;            % vector
%         cNcP  = cP*Nv';         % scalar
%         PN   = P*Nv';          % scalar
%         NcP  = c*(P*Nv'); % scalar - avoid rounding error
%         KNcP = K*NcP;          % scalar
        cKPN = c*K*(P*Nv');          % scalar position . opticdirection

        Ao = -c*(n2 + K*RN.^2);
%         Bo = rN - c*sum(P.*direction,2) - cKNP.*rN;
        Bo = RN .* (1 - cKPN) - c*sum(P.*direction,2);
        Co = c*sum(P.*P,2) + (P * Nv').*(cKPN-2);
    end

    Wo = Bo.^2+Ao.*Co;

    % check if RSI exist
    RSI = (Wo >= 0);
    valid = valid & RSI;

    Wo = sqrt(Wo);
    sw = sign(Bo);
    sw(sw == 0)=1;
    Qo = Bo + sw.*Wo;

    decision = xor(Bo>=0,Outside);

    %suppress divide-by-zero

    if isSphere %Ao is a scalar
        t1 = valid & decision;
        t1(t1) = abs(Ao) > optPrecision;

        t2 = valid & ~decision;
        Qot2 = Qo(t2);
        t2(t2) = abs(Qot2) > optPrecision;

        L=zeros(nRays,1); % gpuArray.
        L(t1) = -Qo(t1)./Ao;
        L(t2) = Co(t2)./Qot2;   
    else %Ao varies across the field
        t1 = valid & decision;
        Aot1 = Ao(t1);
        t1(t1) = abs(Aot1) > optPrecision;

        t2 = valid & ~decision;
        Qot2 = Qo(t2);
        t2(t2) = abs(Qot2) > optPrecision;

        L=zeros(nRays,1); %gpuArray.
        L(t1) = -Qo(t1)./Aot1;
        L(t2) = Co(t2)./Qot2;   
    end

    newPosition = position + bsxfun(@times,direction,L);
    P = bsxfun(@minus,newPosition,surface.position);
%     cP = c*P;
% c is usually small, so multiplying it by P early removes precision in the
% ray position. P is multiplied by a unitary vector (Nv), which also can
% have 1 or 2 elements that are close to zero. Multiplying Nv*P first, 
% and then scaling by c, preserves more precision.el
    if isFlat && ~isNotSimple % [1x3], same for every ray
        N = Nv;
    elseif isSphere % [M x 3]
        N = bsxfun(@minus,Nv,c*P);
    else % [M x 3]
        N = bsxfun(@times,Nv,(1-c*K*(Nv*P')'))-c*P;
    end

    %asphere terms or deformation; iterative solution for both elevation
    %and the surface angle. You could do asphere=[x] to see if it matches
    %a parabola exactly.
    if isNotSimple 
        oldP = P;
        c2 = c*c;
        k1c2=(K + 1)*c2;
        precision2 = optPrecision^2;
        for j=1:maxIterations
            elevation = (Nv * P')';
            % tangentProjection is the projection onto the local
            % coordinate system
            tangentProjection = P - bsxfun(@times,elevation,Nv);
            s2 = sum(tangentProjection.^2,2); % dot(projection,projection)
            if j==1
                z = elevation;
            else % base elevation from tangent plane 
                if c == 0 % plane
                    z = zeros(size(elevation)); % 0 * elevation;
                else
                    z = c*s2./(1+sqrt(1-(k1c2*s2)));
                end
                % P is the nearest point on the parent conic
                P = tangentProjection + z * Nv;
%                 cP = c*P; %local curvature vector, handy in calculating a normal
                N = bsxfun(@times,Nv,(1-c*K*(Nv*P')'))-c*P;
                % N is the local surface-normal vector, not unitary
            end            

            sn = sgn; %concave or convex? establishes sign of s^n
                      %s2 is the square of the radius (aspheres are
                      %axially symmetric)
            Pnear = P;
            % iterative corrections: change Pnear and N
            % the height change is dz. The radial derivative is used to
            % adjust the normal vector
            
            if isASphere
                %aspheric component of curvature, for adding to the normal
                %vector  z(s) = a2 s^2 + a4 s^4 + a6 s^6 + ...
                %    dz/ds(s) = s * (2*a2 + 4*a4*s^2 + 6*a6*s^4)
                %         (one factor of s missing for now. Save it for 
                %         later, and it provides direction
                f = 2;
                dz = zeros(nRays,1); 
                radialDerivative = zeros(nRays,1); 
                for a = surface.asphere
                  radialDerivative = radialDerivative + f * a * sn;
                  sn = sn .* s2;
                  dz = dz + a * sn;
                  f = f + 2;
                end

                % Nearest point on the asphere to the last estimate
                Pnear = tangentProjection + bsxfun(@times,(dz+z),Nv);
                N = N - bsxfun(@times,(Nv * N')' .* radialDerivative,tangentProjection);
               % update the normal vector. 
                % add the radial aspheric component of the derivative to N
                % scale the aspheric component to match the conic normal magnitude
                % This may appear to be bonkers. For a surface defined by
                % the locus of points P(x,y) = [x, y, z(x,y)], the surface
                % normal direction is always [dz/dx, dz/dy, 1]. You can
                % derive it from N = cross(dP/dx, dP/dy). With linear
                % perturbations, P(x,y) = P_conic(x,y) + P_asphere(x,y) we
                % get N(x,y) = N_conic(x,y) + N_asphere(x,y). That's why
                % it's important not to normalize N until all the
                % perturbations are calculated.
            end
            
            % deformations relative to the aperture, not the parent surface
            if isZernike || isDistorted
                if isfield(surface,'center')
                    uv = tangentProjection * surface.local' - surface.center(1:2);     
                elseif isfield(surface,'local')
                    uv = tangentProjection * surface.local';
                else
                    error('a surface needs a local coordinate system to evaluate Zernike coefficients');
                end
 
            if isZernike % sym
               c = surface.zernike;
                if isnumeric(surface.aperture)
                    ap = max(surface.aperture); %circle or annulus 
                else
                    ap = surface.aperture.radius(end); %new form- structure describes aperture. much better
                end
                [dz dxdy] = zernikeXY(c, uv, ap);
                Pnear = Pnear + dz .* Nv;
                % surface.local: [2x3]
                % N: [Nx3]
                % dxdy: [Nx2}      
                N = N - dxdy * surface.local;            
            end
                 
            if isDistorted
                % interpolate a triangulated mesh, where Pnear moves
                % with the interpolated average of the distortion
                % points, and N moves with the interpolated difference

                % vertex indices and barycenters, in local coordinates
%              uv =    tangentLocals = tangentProjection * surface.local';

                % have to iterate over the mesh because a point can be
                % in more than one triangle if it's at a vertex
                for i=1:size(N,1)
                    [triangle,barycenter] = pointLocation(surface.deformation.triangulation,uv(i,:));
                    triVals = surface.deformation.displacement(surface.deformation.triangulation(triangle,:),:);
                    dP = barycenter * triVals; %(bc',triVals')';
                    % dP = mean(dP,1); %multiple triangles get averaged
                    % translates the surface intersection point. 
                    % pointLocation used to return more than one triangle

                    % if deformations are surface-normal, instead of
                    % tangent-normal, you would do dP * normr(N) here
                    Pnear(i,:) = Pnear(i,:) + dP * Nv;

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
                 
% there's a pathology that if a point falls on an edge or vertex, the path length is
% correct, but the tilt comes from only one of the triangles associated
% with that vertex. If any of the barycenters are less than a precision
% threshold, average the tilts from the adjacent triangles as well.
                    notIt = barycenter < optPrecision;
                    s = sum(notIt);
                    if s == 2 %vertex - single
                        vertex = surface.deformation.triangulation.ConnectivityList(triangle,~notIt);
                        dxdy = surface.deformation.vertexSlope(vertex,:);
                    elseif s == 1 %edge - 2 vertices are close
                        vertex = surface.deformation.triangulation.ConnectivityList(triangle,~notIt);
                        edges = all(surface.deformation.edges == sort(vertex),2);
                        %weighted average of vertices? No, better to use
                        %the average of the adjacent triangles
                        %dxdy = sum(surface.deformation.edgeSlope(vertex,:).*barycenter(~notIt),1); 
                        dxdy = surface.deformation.edgeSlope(edges,:);
                    else % look up the slope for that particular vertex
                        dxdy = surface.deformation.slope(triangle,:);
                    end
                    N(i,:) = N(i,:) - dxdy * surface.local;
                end
                end
            end
            end
            % projection and tangentProjection are separate to make
            % it easier to validate while debugging
            projection = tangentProjection;
            % Now project the original ray onto that plane and measure
            % the distance from the intersection point to Pnew
            dL = sum((Pnear - oldP).* N,2) ./ sum(direction.*N,2);

            P = oldP + bsxfun(@times,dL,direction);

            errest = P - Pnear; % distance from the current guess along the plane
            errest = sum(errest.^2,2);
            if all(errest(valid) < precision2) % errest should be identically 0 for on-axis rays
               break
            end
        end
        newPosition = bsxfun(@plus,P,surface.position);
        L = L + dL;%*n2; % change L by the difference between the conic and the distorted surface,
        % remembering that L is in units of direction, which is scaled
        % by the distance vector. we scale dL by n2 later.
    else % not an asphere. Really simple.
        elevation = (Nv * P')';
        projection = P - bsxfun(@times,elevation,Nv);
        if isfield(surface,'local')
            projection = projection * surface.local';
        end
    end

    N = N./sqrt(sum(N.^2,2));
%     N = normr(N);
end

%%
if isfield(surface,'cuy') && surface.cuy < 0
    elevation = -elevation;
end

%% Map onto pixels
if isfield(surface,'pixelSize')
    % [x1 x2 x2 ... ; y1 y2 y3 ...]
%     rays.pixels = round(fix(2*(surface.local * projection')/surface.pixelSize)/2)';
    rays.pixels = round(fix(2*(projection)/surface.pixelSize)/2)';
    rays.pixelSize = surface.pixelSize;
end

%% Aperture test, in tangent plane space. Not optimized for
% computational efficiency because there are some extra
% translations into and out of the local coordinate frame
if optAperture && isfield(surface,'aperture') && ~isempty(surface.aperture) ||isfield(surface,'diameter')
    if isfield(surface,'diameter') && ~isfield(surface,'aperture')
        surface.aperture = surface.diameter / 2;
    end
    % the segment center is defined to 
    if isstruct(surface.aperture)
         if isfield(surface.aperture,'origin')
%              uv = bsxfun(@minus,projection * surface.local',surface.aperture.origin(1:2));
             uv = bsxfun(@minus,projection,surface.aperture.origin(1:2));
         elseif isfield(surface,'center')
%             uv = bsxfun(@minus,projection * surface.local',surface.center(1:2));       
            uv = bsxfun(@minus,projection,surface.center(1:2));       
         elseif isfield(surface,'local')
%             uv = projection * surface.local';
            uv = projection;
         else
            uv = projection;
         end
 
         if strcmp(surface.aperture.type,'pie') % gap, edges, origin, radius
            % [N x 3] * [3 x 2] - [1 x 2]
            % pie slice is with respect to a center and 
            % winding order is u*y - v*x. [y1 y2; -x1 -x2] = [0 1;-1 0] * [x1 x1; y2 y2]' 
            r2 = sum(uv.^2,2); % circular aperture; use winding order, same side of origin
            edges = normr(surface.aperture.edges)'; %transpose 
            e = [0 1;-1 0]*edges;
            g = [1 0;0 -1]*surface.aperture.gap(:);
            % direction is dot product
            side = sign(bsxfun(@plus,uv * e,g')); %which side of the edges? make g a column vector 
            dir = sign(uv * edges) >= 0; % in the right quadrant?
            % side could both be zero (on the line) at a vertex.
            inAP = ((side(:,1) ~= side(:,2))|(side(:,1)==0)) & dir(:,1) & dir(:,2) & ...
                surface.aperture.radius(1)^2 <= r2 & ...
                r2 <= surface.aperture.radius(2)^2; % outer radius
            if isfield(surface.aperture,'parentRadius')
                %w.r.t. segment shape, not parent aperture
                uv = bsxfun(@minus,projection * surface.local',surface.center(1:2));       
                inAP = and(inAP, sum(uv.^2,2) < surface.aperture.parentRadius.^2);
            end
        elseif strcmp(surface.aperture.type,'circle')
            r2 = sum(uv.^2,2); % circular aperture
            inAP = r2 <= surface.aperture.radius^2;
        elseif strcmp(surface.aperture.type,'annulus')
            r2 = sum(uv.^2,2); % circular aperture
            inAP = surface.aperture.radius(1)^2 <= r2 & r2 <= surface.aperture.radius(2)^2;          
        elseif strcmp(surface.aperture.type,'rectangle')
            % surface.aperture.bounds = |u1 u2| 
            %                           |v1 v2|
            inAP = surface.aperture.bounds(1) <= uv(:,1) & uv(:,1) <= surface.aperture.bounds(3)& ...
                surface.aperture.bounds(2) <= uv(:,2) & uv(:,2) <= surface.aperture.bounds(4);
        elseif strcmp(surface.aperture.type,'polygon')
            x = surface.aperture.bounds(:,1);
            y = surface.aperture.bounds(:,2);
            if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
                x(end+1)=x(1);
                y(end+1)=y(1);
            end
            [in,on] = inpolygon(uv(:,1),uv(:,2),x,y);
            inAP = in|on;
       else
            error(['invalid aperture type ' surface.aperture.type]);
        end
    else % center, local coordinate sys, radius|rectangle|ellipse|polygon
        if isfield(surface,'center')
            uv = bsxfun(@minus,projection * surface.local',surface.center(1:2));       
        elseif isfield(surface,'local')
            uv = projection * surface.local';
        else
            uv = projection;
        end
        
        M = numel(surface.aperture);
        if M == 4 % rectangle
            % |u1 v1|  contains [u v]
            % |u2 v2|
            inAP = surface.aperture(1) <= uv(:,1) & uv(:,1) <= surface.aperture(2)& ...
                surface.aperture(3) <= uv(:,2) & uv(:,2) <= surface.aperture(4);
        elseif M == 3 % ellipse
            error('Invalid aperture shape')
        elseif M < 3
            r2 = sum(uv.^2,2); % circular aperture
            if M == 2 % [inner outer]
                inAP = surface.aperture(1)^2 <= r2 & r2 <= surface.aperture(2)^2;
            elseif M == 1 % outer
                inAP = r2 <= surface.aperture^2;
            else
                error([label ' aperture must have 1, 2, or 4 elements']);
            end
        else % polygon
            x = surface.aperture(:,1);
            y = surface.aperture(:,2);
            if x(1) ~= x(end) || y(1) ~= y(end) %close the polygon
                x(end+1)=x(1);
                y(end+1)=y(1);
            end
            [in,on] = inpolygon(uv(:,1),uv(:,2),x,y);
            inAP = in|on;
        end
    end
    RSI = RSI & inAP;
end
% No RSI
noRSI = oldValid & ~RSI;
rays.status(noRSI) = -1; % No RSI
valid = valid & RSI;

%% Neg. Ray test
% At this point, valid = oldValid & RSI 
negRays = (L<0);    
rays.status(valid & negRays) = -2; % Negative ray
if ~optNegRays
    valid = valid & ~negRays;
end

%% Found the RSI! Now what? Update the path length
position = newPosition;
distance = n2 .* L;
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
    elseif isfield(surface,'n') 
        ne2 = surface.n*surface.n; 
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
            rays.status(errTIR) = -1;
            valid(errTIR) = false;

%             t = bsxfun(@times,N,np-surface.*sqrt(TIR));
            t = N .* (np-sgn.*sqrt(TIR));
            direction(valid,:) = direction(valid,:) - t(valid,:);

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


if isfield(surface,'display'), rays.display = surface.display; end    
if isfield(rays,'opl'), rays.opl = rays.opl + distance; else rays.opl = distance; end

% 
%todo: transform the coordinate system that the rays drag along with them
% u = new direction
% v = original direction
% If u . v < 1, they're in opposite directions, so need a sign
if isfield(rays,'local') 
    local = rays.local';
    if size(local,2) == 1
        local = repmat(local,1,nRays);
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
                local(:,i) = -(r*(local(:,i)));
            else
                local(:,i) = (r*(local(:,i)));
            end
        end
    end
    rays.local = local';
end

rays.position = position;
rays.direction = direction;
rays.status = status;
rays.valid = valid;
rays.normal = N;
rays.surface = surface;
if optDebug
    projection(~valid,:)=0;
    rays.projection = projection; %in the global coordinate system, relative to the vertex position
    rays.distance = distance;
    rays.elevation = elevation;
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