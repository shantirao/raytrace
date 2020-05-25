function glass = makeGlass(catalog,name,wavelengthUnits)
% catalog is schoot, ohara, or hoya
% name is the glass name
% wavelengthUnits defaults to mm
glass = struct();

if nargin < 3
    wavelengthUnits = 'mm';
end

if strcmp(catalog,'schott')    %Sellmeier equation
    data = lookup(importTextData('Schott_2013a.txt'),name);
    glass.B1=data{1};
    glass.B2=data{2};
    glass.B3=data{3};
    glass.C1=data{4};
    glass.C2=data{5};
    glass.C3=data{6};
    glass.units='um';
    glass.evaluate=@glassSellmeier;
elseif strcmp(catalog,'ohara')    %Sellmeier equation
    data1 = lookup(importTextData('OHARA_2012A.txt'),name);
    data2 = lookup(importTextData('OHARA_2012B.txt'),name);
    glass.B1=[data1{1},data2{1}];
    glass.B2=[data1{2},data2{2}];
    glass.B3=[data1{3},data2{3}];
    glass.C1=[data1{4},data2{4}];
    glass.C2=[data1{5},data2{5}];
    glass.C3=[data1{6},data2{6}];
    glass.units='um';
    glass.range=.001129;
    glass.evaluate=@glassSellmeier;
elseif strcmp(catalog,'hoya')    %Schott equaiton
    data = lookup(importTextData('HOYA_2013.txt'),name);
    glass.A0=data{1};
    glass.A1=data(2);
    glass.A2=data{3};
    glass.A3=data{4};
    glass.A4=data{5};
    glass.A5=data{6};
    glass.units='um';
    glass.evaluate=@glassSchott;
end

function coeffs = lookup(data,name)
for i=1:numel(data)
    if strcmp(name,data{i}{1})
        coeffs = {data{i}{2:end}};
    end
end

function [data, header] = importTextData(filename)

%% Open the text file.
fid = fopen(filename,'r');

data = {};
header={};
line = fgetl(fid);
while ischar(line)
    if line(1) ~= '&' && line(1) ~= '#' && line(1) ~= '%' && numel(line)>0
        data{end+1} = textscan(line,'%s %f %f %f %f %f %f');
    else
        header{end+1}=line;
    end
    line = fgetl(fid);
end
fclose(fid);

