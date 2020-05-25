function index2 = glassSchott(glass,wavelength, units)
% glass is an object with fields B1, B2, B3
% evaluates the Sellmeier emperical dispersion formula
% takes wavelength in microns, returns index^2

if nargin > 2
    scale = scaleUnits(glass.units)/scaleUnits(units);
    wavelength = wavelength * scale;
end

L2 = wavelength*wavelength;
L4 = L2 * L2;
L6 = L4 * L2;
L8 = L4 * L4;
index2 = glass.A0 + glass.A1*L2 + glass.A2/L2 +  glass.A3/L4 + glass.A4/L6 + glass.A5/L8;

end