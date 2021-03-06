%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/rbirger/Dropbox/Yale/Typhoid Pitzer Lab/AMR/ARVacProject/MatlabCode/Gavi_cov_unc.csv
%
% Auto-generated by MATLAB on 21-Jun-2020 15:29:47

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 12);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Vaccine", "Version", "ISO3", "Gavi73", "Year", "Delivery", "IntroductionDate", "GaviSupport", "Population", "FVP", "FVPType", "Coverage"];
opts.VariableTypes = ["categorical", "double", "categorical", "categorical", "double", "categorical", "double", "string", "double", "double", "categorical", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "GaviSupport", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Vaccine", "ISO3", "Gavi73", "Delivery", "GaviSupport", "FVPType"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Version", "TrimNonNumeric", true);
opts = setvaropts(opts, "Version", "ThousandsSeparator", ",");

% Import the data
Gavicovunc = readtable("DataFiles/Gavi_cov_unc.csv", opts);


%% Clear temporary variables
clear opts