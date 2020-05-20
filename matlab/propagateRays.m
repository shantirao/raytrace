function [rays, distance, projection, elevation] = propagateRays(rays, surface, options)
% see raytrace.m for instructions

optAperture = true;
optNegRays = false;
optSegments = true;
optParallel = false;
optPrecision = 1e-14;
optDebug = false;

if nargin > 2
    if isfield(options,'aperture'), optAperture = options.aperture; end
    if isfield(options,'negative'), optNegRays = options.negative; end
    if isfield(options,'segments'), optSegments = options.segments; end
    if isfield(options,'parallel'), optParallel = options.parallel; end
    if isfield(options,'precision'), optPrecision = options.precision; end
    if isfield(options,'debug'), optDebug = options.debug; end
end

maxIterations = 10;
nRays = size(rays.position,1);
position = (rays.position); % could be a gpuArray
% note that rays.direction is scaled by index^2 
direction = (rays.direction); % could be a gpuArray
if ~isfield(rays,'status'), rays.status = ones(nRays,1); end
if ~isfield(rays,'valid'), rays.valid = true(nRays,1); end
if isfield(rays,'n2'), n2 = rays.n2; else n2=1; end
if isfield(rays,'n'), n = rays.n; else n=1; end
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
isFlat = ~isfield(surface,'cuy') || abs(surface.cuy) < optPrecision;
isSphere = ~isFlat && (~isfield(surface,'K') || abs(surface.K) < optPrecision);
isASphere = isfield(surface,'asphere') ;
isDistorted = isfield(surface,'deformation');
isNotSimple = isASphere || isDistorted;

%in this branch, N is the surface normal vector. For planes it's a
%single vector. For curved surfaces, it's Nx3, a 3-vector for each
%intersection point.
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

else %not simple: distorted flat, or a conic  
    % Initialize sign conventions
    if isFlat
        c = 0;
        K = 0;
        Outside = false; 
        Nv = surface.direction;
        sgn = -1;
        % Outside means the normal vector points to the outside
    else
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
        else
            Nv = surface.direction;
            c  = surface.cuy;
            sgn = -1;
        end            
    end

    P = bsxfun(@minus,position,surface.position);

    if isSphere % K = 0 simplifies the equations
        Ao = -c*n2;
        Bo = (direction*Nv') - sum((c*P).*direction,2);
        Co = c*sum(P.^2,2) - 2*(P * Nv');
    else % K = anything; reduces to the previous case for K = 0 or a flat for cuy = 0
        rN   = direction*Nv';
        cP   = c*P;
        NcP  = cP*Nv';
        NcPK = K*NcP;

        Ao = -c*(n2 + K*rN.^2);
        Bo = rN - sum(cP.*direction,2) - NcPK.*rN;
        Co = sum(cP.*P,2) + (P * Nv').*(NcPK-2);
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
    cP = c*P;
    if isFlat && ~isNotSimple % [1x3], same for every ray
        N = Nv;
    elseif isSphere % [M x 3]
        N = bsxfun(@minus,Nv,cP);
    else % [M x 3]
        N = bsxfun(@times,Nv,(1-K*(Nv*cP')'))-cP;
    end

    %asphere terms or deformation; iterative solution
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
            s2 = sum(tangentProjection.^2,2);
            if j==1
                z = elevation;
            else
                if c == 0 % plane
                    z = 0 * elevation;
                else
                    z = c*s2./(1+sqrt(1-(k1c2*s2)));
                end
                % P is the nearest point on the parent conic
                P = tangentProjection + z * Nv;
                cP = c*P; %local curvature vector, handy in calculating a normal
                N = bsxfun(@times,Nv,(1-K*(Nv*cP')'))-cP;
                % N is the local surface-normal vector
            end            

            sn = sgn;
            radialDerivative = zeros(nRays,1); 
            Pnear = P;
            %iterative corrections: change Pnear and N
            if isASphere
                %aspheric component of curvature, for adding to the normal
                %vector  z(s) = a2 s^2 + a4 s^4 + a6 s^6 + ...
                %       dz(s) = 2*a2 + 4*a4*s^2 + 6*a6*s^4
                %this would be a good place for Zernike perturbations
                f = 2;
                dz = 0;
                for a = surface.asphere
                  radialDerivative = radialDerivative + f * a * sn;
                  sn = sn .* s2;
                  dz = dz + a * sn;
                  f = f + 2;
                end

                % Nearest point on the asphere to the last estimate
                Pnear = tangentProjection + bsxfun(@times,(dz+z),Nv);
                % update the normal vector. There is a 
                % add the radial aspheric component of the derivative to N
                % scale the aspheric component to match the conic normal magnitude
                N = N - bsxfun(@times,(Nv * N')' .* radialDerivative,tangentProjection);
            end

            if isDistorted
                % interpolate a triangulated mesh, where Pnear moves
                % with the interpolated average of the distortion
                % points, and N moves with the interpolated difference

                % vertex indices and barycenters, in local coordinates
                tangentLocals = tangentProjection * surface.local';

                % have to iterate over the mesh because a point can be
                % in more than one triangle if it's at a vertex
                for i=1:size(N,1)
                    [vi,bc] = pointLocation(surface.deformation.triangulation,tangentLocals(i,:));
                    triVals = surface.deformation.displacement(surface.deformation.triangulation(vi,:),:);
                    dP = bc * triVals; %(bc',triVals')';
                    dP = mean(dP,1); %multiple triangles get averaged
                    % translates the surface intersection point. 
                    Pnear(i,:) = Pnear(i,:) + dP;

                    % the surface normal moves a bit, too. It can be
                    % expressed as a rotation. Use the curl of the
                    % triangular facets to rotate the normal vector
                    % if more than one triangle, take the average curl

                    rotvec = mean(surface.deformation.curl(vi,:),1);
                    Q = rotationMatrix(normr(rotvec),norm(rotvec));
                    N(i,:) = Q * N(i,:)';
                end
            end

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
    end

    N = normr(N);
end

% pixels
if isfield(surface,'pixelSize')
    % [x1 x2 x2 ... ; y1 y2 y3 ...]
    rays.pixels = round(fix(2*(surface.local * projection')/surface.pixelSize)/2)';
end

% Aperture test, in tangent plane space. Not optimized for
% computational efficiency because there are some extra
% translations into and out of the local coordinate frame
if optAperture && isfield(surface,'aperture') && ~isempty(surface.aperture)        
    % the segment center is defined to 
    if isstruct(surface.aperture)
        if strcmp(surface.aperture.type,'pie')
            uv = bsxfun(@minus,projection * surface.local',surface.aperture.origin(1:2)); 
            % pie slice is with respect to a center and 
            r2 = sum(uv.^2,2); % circular aperture; use winding order, same side of origin
            % winding order is u*y - v*x. [y1 y2; -x1 -x2] = [0 1;-1 0] * [x1 x1; y2 y2]' 
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

% Neg. Ray Length
% At this point, valid = oldValid & RSI 
negRays = (L<0);    
rays.status(valid & negRays) = -2; % Negative ray
if ~optNegRays
    valid = valid & ~negRays;
end

position = newPosition;
distance = n2 * L;
%     opl = opl + distance;

surfType = 0; % stop or reference
if isfield(surface,'type') 
    surfType = surface.type;
end
%reflect or refract? Note N is normalized
switch surfType
    case 1 %reflect
        if isFlat && ~isNotSimple
            direction = bsxfun(@minus,direction,2*(N*direction')'*N);
        else % dot product = sum(N.*direction,2)
            delta = 2*bsxfun(@times,N,sum(N.*direction,2));
            direction = bsxfun(@minus,direction,delta);
        end        
    case 2 %refract
        ne2 = surface.n*surface.n;
        np=sum(bsxfun(@times,N,direction),2);
        surface = sign(np);
        surface(surface==0) = 1;
        TIR = ne2 - n2 + np.^2;

        %error detection
        errTIR = TIR<0;
        lost = valid;
        lost(valid) = errTIR;
        rays.status(lost) = -1;
        valid(lost) = false;

        t = bsxfun(@times,N,np-surface.*sqrt(TIR));
        direction(valid,:) = direction(valid,:) - t(valid,:);

        % new surface
        rays.n2 = ne2;
%         rays.n = surface.n;

    case 4 %polarization: reflect and possibly refract
        % make a new polarizationReflection() function
        if isfield(rays,'local') %project the TE vector onto the reflection plane
            TE = rays.local;
        else
            TE = bsxfun(@cross,dirnorm,N);
        end
        dirnorm = direction / n;
        s = bsxfun(@cross,dirnorm,N); %ray direction before reflection [N x 3]
        s = normr(s); % problem at normal incidence
        cosphi = sum(s.*rays.local,2)/n; % dot product, scaled down by the incident ray 
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

if isfield(surface,'display'), rays.display = surface.display; end    
if isfield(rays,'opl'), rays.opl = rays.opl + distance; else rays.opl = distance; end

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