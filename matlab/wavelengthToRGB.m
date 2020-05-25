function rgb = wavelengthToRGB(wavelength, gamma)
%http://www.noah.org/wiki/Wavelength_to_RGB_in_Python

if nargin < 2
    gamma = 0.8;
end

rgb = [];
for i=1:numel(wavelength)
    %convert from mm to nm;
     lambda = wavelength(i) * 1e6;
    if lambda >= 380 && lambda <= 440
        attenuation = 0.3 + 0.7 * (lambda - 380) / (440 - 380);
        R = ((-(lambda - 440) / (440 - 380)) * attenuation) .^ gamma;
        G = 0.0;
        B = (1.0 * attenuation) .^ gamma;
    elseif lambda >= 440 && lambda <= 490
        R = 0.0;
        G = ((lambda - 440) / (490 - 440)) .^ gamma;
        B = 1.0;
    elseif  lambda >= 490 && lambda <= 510
        R = 0.0;
        G = 1.0;
        B = (-(lambda - 510) / (510 - 490)) .^  gamma;
    elseif  lambda >= 510 && lambda <= 580
        R = ((lambda - 510) / (580 - 510)) .^  gamma;
        G = 1.0;
        B = 0.0;
    elseif  lambda >= 580 && lambda <= 645
        R = 1.0;
        G = (-(lambda - 645) / (645 - 580)) .^  gamma;
        B = 0.0;
    elseif  lambda >= 645 && lambda <= 750
        attenuation = 0.3 + 0.7 * (750 - lambda) / (750 - 645);
        R = (1.0 * attenuation).^  gamma;
        G = 0.0;
        B = 0.0;
    else
        R=0;
        G=0;
        B=0;
    end

    rgb(end+1,:) = [R,G,B];
end

end
