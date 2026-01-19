% Import numerical HHG output data from a text file.
function photo = importfile_hhg_data(filename, dataLines)

if nargin < 2
    dataLines = [1, Inf];
end


opts = delimitedTextImportOptions("NumVariables", 1);

opts.DataLines = dataLines;
opts.Delimiter = ",";

opts.VariableNames = "VarName1";
opts.VariableTypes = "double";

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

photo = readtable(filename, opts);

photo = table2array(photo);
end