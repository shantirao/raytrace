function [start, finish] = importVV(filename, startRow, endRow)
%importVV Import numeric data from a text file as a matrix.
%   [startRays,finishRays] = loadExample(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   EXAMPLE001VV = loadExample(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%   [ X(0), Y(0), Z(0) ] , [ L(0) , M(0) , N(0) ] , [ X(n), Y(n), Z(n) ] , [ L(n) , M(n) , N(n) ] , OPL
% Example:
%   example001VV = importfile('example_001_VV.txt', 3, 3212);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2014/03/03 11:20:55

%% Initialize variables.
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%24f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
ex = [dataArray{1:end-1}];

start.N = size(ex,1);
start.chief=1;
start.position = ex(:,1:3);
start.direction = ex(:,4:6);
start.opl = zeros(start.N,1);

finish.N = start.N;
finish.chief=1;
finish.position = ex(:,7:9);
finish.direction = ex(:,10:12);
finish.opl = ex(:,13);