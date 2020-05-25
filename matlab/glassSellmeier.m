function index2 = glassSellmeier(glass,wavelength, units)
% glass is an object with fields B1, B2, B3
% evaluates the Sellmeier emperical dispersion formula
% takes wavelength in microns, returns index^2

if nargin > 2
    scale = scaleUnits(glass.units)/scaleUnits(units);
    wavelength = wavelength * scale;
end

index2 = 1 + glass.B1(1)*wavelength.^2./(wavelength.^2-glass.C1(1))+glass.B2(1)*wavelength.^2./(wavelength.^2-glass.C2(1))+glass.B3(1)*wavelength.^2./(wavelength.^2-glass.C3(1));

if isfield(glass,'range')
    a = 1 + glass.B1(2)*wavelength.^2./(wavelength.^2-glass.C1(2))+glass.B2(2)*wavelength.^2./(wavelength.^2-glass.C2(2))+glass.B3(2)*wavelength.^2./(wavelength.^2-glass.C3(2));
	index2(wavelength > glass.range) = a(wavelength > glass.range);
end

end