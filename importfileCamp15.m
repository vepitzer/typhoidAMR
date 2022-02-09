function Camp15Impact = importfileCamp15(filename, dataLines)
%IMPORTFILE Import data from a text file
%  CAMP15IMPACT73TOTAL = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  CAMP15IMPACT73TOTAL = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Camp15Impact73Total = importfile("Camp15Impact_73_Total.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 13-Aug-2020 19:32:14

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["country", "CasesAvertedMedian", "CasesAvertedMin", "CasesAvertedMax", "DeathsAvertedMedian", "DeathsAvertedMin", "DeathsAvertedMax", "DALYsAvertedMedian", "DALYsAvertedMin", "DALYsAvertedMax"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "country", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "country", "EmptyFieldRule", "auto");

% Import the data
Camp15Impact = readtable(filename, opts);

end