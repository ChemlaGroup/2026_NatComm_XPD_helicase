% 230209, AVT

% Load all datasets of interest into a single cell array, in one fell
% swoop. 
% In resulting 'data' cell array: data{x}{y} corresponds to trace number
% 'y' for the construct indexed by 'x' (with construct indeces as defined
% in config file). Each such element data{x}{y} is a structure containing
% various kinds of data associated with the trace.

% Sample usage:
% data = LoadData(0);
% data = LoadData(0, [9 14 15 16 17 24 19 20 21 22 23 27]);
% dataq = LoadData(1);

function data= LoadData(baselineselect, constructs)

if nargin<2
    %constructs = [7 8 9 10 11];
    constructs = [9 14 15 16 17 18 19 20 21 22 23];
    %constructs = [19 20 21 23];
end

data={}; 

for c=constructs
  [data{c}, ~] = DatasetLoading_Baselineselect(c, baselineselect);
end

end