% -------------------------------------------------------------------------------
% [psf] = MakeBroadbandPSF(amp, opd, bandpass, sampling, dfs_tilt, [geom_pupil]);
% -------------------------------------------------------------------------------
%
% Generates a broadband PSF for the specified exit pupil model and bandpass.
%
%
% INPUTS:
%
%      amp - exit pupil amplitude function
%
%      opd - exit pupil optical path difference function
%           (same units as wavelenth and pix_size)
%
% bandpass - spectral range and relative response of detection process
%
%            (e.g.)       Wavelength           Intensity Xmission
%                   ----------------------------------------------
%                   bandpass (1,1) = 631e-9    bandpass(1,2) = 0.5
%                   bandpass (2,1) = 633e-9    bandpass(2,2) = 1.0
%                   bandpass (3,1) = 635e-9    bandpass(3,2) = 0.5
%
% sampling - sampling parameters = [ Nyquist Fraction   [Wavelength]]
%
%            (e.g.) sampling = [2]      -> 2x Nyquist at bandpass centroid
%                   sampling = [2 1e-6] -> 2x Nyquist at wavelength=1e-6
%
% dfs_tilt = [kx [ky]]  - used to simulate DFS images using fixed tilt
%
%            (e.g.) dfs_tilt = 0                ->  Normal PSF generation
%                   dfs_tilt = 20e-6 * [1 .5]   ->  DFS on a [1 .5] orientation
%
% geom_pupil = (optional) geometric pupil mask which defines the pupil extent
%              to setup sampling.  If ommitted the default is to establist the
%              sampling based upon the non-zero extent of pupil amplitude
%              function
%
% OUTPUTS:
%
%      psf = The broadband PSF
% bandpass = The true bandpass represented in the PSF
%
% -------------------------------------------------------------------------------
% Version Nov-30-2004:  Joseph Green, Jet Propulsion Laboratory
% -------------------------------------------------------------------------------

function [psf, bandpass] = MakeBroadbandPSF(amp, opd, bandpass, sampling, dfs_tilt, geom_pupil)


    % ------------------------------------------------------------------
    % Pad Optical bandpass by a zero response on either side
    % ------------------------------------------------------------------
    bandpass = sortrows(bandpass);

    wl  = bandpass(:,1);
    res = bandpass(:,2);

    if length(wl)>1
       dwl  = (max(wl)-min(wl))/(length(wl)-1);
    else
       dwl  = wl/100;
    end

    wl   = [(wl(1)-dwl); wl ; (wl(length(wl)) + dwl)];
    res  = [0; res; 0];

    % ------------------------------------------------------------------
    % Setup the pupil padding as a function of lambda
    % ------------------------------------------------------------------

    if nargin < 6
       geom_pupil = amp~=0;
    end

    [xmin xmax ymin ymax] = find_boundary (geom_pupil, .5);

    %BK
    npix_pup = max(xmax-xmin,ymax-ymin) + 1;
%    npix_nyq = 2*npix_pup - 1;
    npix_nyq = 2*npix_pup;

    if length(sampling) == 1
       lambda_o = sum (wl.*res)./sum(res);
    else
       lambda_o = sampling(2);
       sampling = sampling(1);
    end

    npix_pad = round(npix_nyq * sampling * wl ./ lambda_o);
    %npix_pad = npix_pad - mod(npix_pad,2);
    %BK
    npix_pad = unique(npix_pad);

    % ------------------------------------------------------------------
    % Recompute the true wavelenths represented by the padding
    % ------------------------------------------------------------------

    lambda   = lambda_o * npix_pad ./ (npix_nyq * sampling);

    % ------------------------------------------------------------------
    % Interpolate the bandpass response over the recomputed wavelengths
    % ------------------------------------------------------------------

    flux = res / dwl;

    if min(lambda) < 0
        error('Wavelengths must be positive');
    end
    l = [];
    r = [];
    p = [];
    for n=2:length(lambda)-1

       dl1 = ( lambda(n)   - lambda(n-1) )/2;
       dl2 = ( lambda(n+1) - lambda(n)   )/2;

       l0 = lambda(n);
       l1 = lambda(n) - dl1;
       l2 = lambda(n) + dl2;

       f0 = interp1(wl,flux,l0);
       if isnan(f0); f0=0; end

       f1 = interp1(wl,flux,l1);
       if isnan(f1); f1=0; end

       f2 = interp1(wl,flux,l2);
       if isnan(f2); f2=0; end

       r1 = (f0+f1)/2 * dl1;
       r2 = (f0+f2)/2 * dl2;

       r0 = r1+r2;

       r(end+1) = r0; % = [r r0];
       l(end+1) = lambda(n); % = [l lambda(n)];
       p(end+1) = npix_pad(n); % = [p npix_pad(n)];

    end

    lambda   = l;
    res_int  = r * sum(res)/sum(r);
    npix_pad = p;


    if 1==0
        res_int    = interp1(wl,res,lambda);
        v          = find(isnan(res_int));
        res_int(v) = 0;
    end

    % ------------------------------------------------------------------
    % Setup the DFS tilts
    % ------------------------------------------------------------------

    if dfs_tilt ~= 0
        z2 = zernike_mode (geom_pupil~=0,2);
        z3 = zernike_mode (geom_pupil~=0,3);

        dfs_x = dfs_tilt(1) * (1 - lambda./mean(lambda));
        dfs_y = 0*dfs_x;
        if max(size(dfs_tilt))>1
           dfs_y = dfs_tilt(2) * (1 - lambda./mean(lambda));
        end
    else
        dfs_x = 0;
        dfs_y = 0;
        z2 = 0;
        z3 = 0;
    end
    % ------------------------------------------------------------------
    % Build Up the the broad-band PSF image
    % ------------------------------------------------------------------

    ri   = find(res_int>0);
    npix = max (npix_pad(ri));
    Nmax = length(ri);
    bandpass = zeros(Nmax,2);

%     disp('MakeBroadbandPSF: Generating PSF...');
% parallelize on the CPU if it's too big for a GPU
    if ~exist("gpuArray","builtin") || gpuDeviceCount == 0 || Nmax > 8 || sampling > 2 % delete(gcp) to close the current parallel pool
        psf = zeros(npix,npix);
        if dfs_tilt ~= 0
            parfor n=1:length(ri)
                k = ri(n); %skips the 0-weighted elements in the bandpass.
    %        disp (sprintf(' %3d/%3d  Nyquist Factor=%5.4f  lambda=%9.7g  response=%9.7g',n,length(ri),npix_pad(k)/npix_nyq,lambda(k),res_int(k)));
               p = amp .* exp(1i*(2*pi/lambda(k)) * (k + z2*dfs_x(k) + z3*dfs_y(k)) );
               p = ifftshift(resizeAndCenter(p, npix_pad(k))); %quadswap(pad(p, npix_pad(k)), -1);
%                p = quadswap(pad(p, npix_pad(k)), -1);
               P = abs(fftshift(fft2(p))).^2;
               P = resizeAndCenter(P * (res_int(k) / sum(P(:))), npix);
               psf(:,:,n) =  P;
               bandpass(n,:) = [lambda(k) res_int(k)];
            end
        else
            parfor n=1:length(ri)
                k = ri(n); %skips the 0-weighted elements in the bandpass.
                p = amp .* exp(1i*(2*pi/lambda(k)) * opd);
                p = ifftshift(resizeAndCenter(p, npix_pad(k))); %quadswap(pad(p, npix_pad(k)), -1);
      %          P = abs((fft2(p))).^2;
                P = abs(fftshift(fft2(p))).^2;
                P = resizeAndCenter(P * (res_int(k) / sum(P(:))), npix);
                psf(:,:,n) =  P;
                bandpass(n,:) = [lambda(k) res_int(k)];
            end
        end
        psf=sum(psf,3); %collapse
    else
%         lambda = gpuArray(lambda);
         opd = gpuArray(opd);
         amp = gpuArray(amp);
         res_int = gpuArray(res_int);
        psf = zeros(npix,npix,'gpuArray');
        if dfs_tilt ~= 0
            dfs_x = gpuArray(dfs_x);
            dfs_y = gpuArray(dfs_y);
            for n=1:length(ri)
                k = ri(n); %skips the 0-weighted elements in the bandpass.
    %        disp (sprintf(' %3d/%3d  Nyquist Factor=%5.4f  lambda=%9.7g  response=%9.7g',n,length(ri),npix_pad(k)/npix_nyq,lambda(k),res_int(k)));
               p = amp .* exp(1i*(2*pi/lambda(k)) * (opd + z2*dfs_x(k) + z3*dfs_y(k)) );
               p = ifftshift(resizeAndCenter(p, npix_pad(k))); %quadswap(pad(p, npix_pad(k)), -1);
               P = abs(fftshift(fft2(p))).^2;
               P = resizeAndCenter(P * (res_int(k)/ sum(P(:))), npix);
               psf = psf + P;
                clear p;
                clear P;
                bandpass(n,:) = [gather(lambda(k)) gather(res_int(k))];
            end
        else
            for n=1:length(ri)
                k = ri(n); %skips the 0-weighted elements in the bandpass.
                p = amp .* exp(1i*(2*pi/lambda(k)) * opd);
                p = resizeAndCenter(ifftshift(p), npix_pad(k)); %p = padarray(p,[padSize padSize]); % pad(p, npix_pad(k))
                P = abs(fftshift(fft2(p))).^2;
%                 P ./ sum(P(:)); % normalize it
 %               padSize = floor((npix/length(P))/2);
%                P = padarray(P,[padSize padSize]);
                P = resizeAndCenter(P * (res_int(k)/ sum(P(:))), npix);
                psf = psf + P;
                clear p;
                clear P;
                bandpass(n,:) = [gather(lambda(k)) gather(res_int(k))];
%                 disp(sprintf('%d %d',n,g.FreeMemory));
            end
        end
        psf = gather(psf);
    end
return;

function [f_out] = quadswap (f_in, direction)

   if direction ~= -1
      f_out = fftshift(f_in);
      return
   end

   f_out = ifftshift(f_in);
   return

   [nrows, ncols]  = size(f_in);
   xc = floor(ncols/2 + 1);
   yc = floor(nrows/2 + 1);
   Q1 = f_in (yc:nrows,xc:ncols);
   Q2 = f_in (yc:nrows,  1:xc-1);
   Q3 = f_in (1:yc-1,    1:xc-1);
   Q4 = f_in (1:yc-1,  xc:ncols);

   f_out = [ [Q1 Q2]; [Q4 Q3] ];

return

% ------------------------------------------------------
% [xmin xmax ymin ymax] = find_boundary (img, tol)
% ------------------------------------------------------
%
% Finds bounding xmin, xmax, ymin and ymax values
% where
%
%    abs(img(ymin:ymax,xmin:xmax)) > tol*max(abs(img))
%
% ------------------------------------------------------

function [xmin, xmax, ymin, ymax] = find_boundary (img, tol)

maxv = max(max(abs(img)));

if maxv == 0
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
    return;
end

mask = (abs(img)./maxv > tol);

[yy, xx] = find (mask == 1);

xmin = min(xx);
xmax = max(xx);
ymin = min(yy);
ymax = max(yy);
return;
