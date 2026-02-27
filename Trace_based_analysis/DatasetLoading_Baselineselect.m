
% Loading the data used in the analysis. Optional: exclude traces failing
% pre-defined noise criteria.

% Needs updating as more data collected and processed, and/or baseline selection criteria changed.

% baselineselect: whether all processed data for construct is returned, or
% the subset passing noise quality cutoffs

% includelist is only non-empty if baselineselect was on

% Sample usage:
% [data7, ~] = DatasetLoading_Baselineselect(7, 0);
% [data14, ~] = DatasetLoading_Baselineselect(14, 1);
% [data9, ~] = DatasetLoading_Baselineselect(9, 1);data{9}=data9;
% [data16M, ~] = DatasetLoading_Baselineselect(160171,0);

% One fell swoop: (or use LoadData)
% constructs = [9 14 15 16]; data={}; dataq={};
% for c=constructs
% [dataq{c}, ~] = DatasetLoading_Baselineselect(c, 1); [data{c}, ~] = DatasetLoading_Baselineselect(c, 0);
% end

function [dataq, includelist] = DatasetLoading_Baselineselect(construct, baselineselect)
% Load data

global analysis
global datadirectories

if nargin < 2
    baselineselect = 1;
end

if nargin <1
    construct = 9; 
end


if construct == 9 % 5' control
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];
    
elseif construct ==14 % 5' abasic
    cd(datadirectories)
    load('sample_traces_abasic_5set'); 
    alldata = [sample_traces_abasic]; 

elseif construct ==15 % 5' on-strand mismatch
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];
    
elseif construct ==16   % 5' CPD
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];
    
elseif construct ==17   % 5' fluorescein
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];
    
elseif construct ==19 % 3' control
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];

elseif construct == 20 % 3' CPD
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];

elseif construct == 21 % 3' Fluorescein
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];
    
elseif construct == 22 % 3' Abasic
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];    
    
elseif construct == 23 % 3' Mismatch   
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];   
    
elseif construct == 24 % 5' Mismatch #2 (off-strand)
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];   

elseif construct == 27 % 3' Mismatch #2 (off-strand) 
    error('No data currently provided for this construct.')
%     cd('DataContainingDirectory/');
%     load('datafile1'); load('datafile2');
%     alldata = [data1 data2];   

end

cd(analysis)

% baseline-select

% Choose subset of total data by amount of baseline noise -- based on third
% baseline, after XPD activity ends.
% Based on code in baseline_dists_calc; see for further comments.

dataq ={};
includelist = []; excludelist = [];
baselinebinw = 0.1;
fit_constraint = 3;

%Baseline quality parameters:
b3_mean_threshold = 2;
b3_stddev_threshold = 1;
b3fit_mean_threshold = 1;
b3fit_stddev_threshold = 0.8; 
%Width-related value (c1) returned by Gaussian fit: actually std. dev*sqrt(2).
% What is called fit_stddev in this code needs to be divided by sqrt(2)
% to get real std. dev.

if baselineselect
    
    for i= 1:length(alldata)
        [b3_mean, b3_stddev, b3fit_mean, b3fit_stddev] = trace_selection_baseline_noise(alldata{i}, baselinebinw, fit_constraint, b3_mean_threshold, b3_stddev_threshold);
        
        if (abs(b3_mean)<=b3_mean_threshold && b3_stddev<=b3_stddev_threshold)
            includelist= [includelist i];
            dataq = [dataq alldata{i}];
        elseif (abs(b3fit_mean) <= b3fit_mean_threshold && b3fit_stddev <=b3fit_stddev_threshold)
            includelist= [includelist i];
            dataq = [dataq alldata{i}];
        else
            excludelist = [excludelist i];
        end
    end
    clear b3_mean b3_stddev b3fit_mean b3fit_stddev;
    Ntraces = length(includelist);
    
    disp(['Number of traces passing cutoff: ' num2str(Ntraces) ' out of ' num2str(Ntraces+length(excludelist))]);
    %disp('Includelist:')
    %includelist
    
else
    disp('NO BASELINE SELECTION!')
    dataq = alldata;
    Ntraces = size(dataq, 2);
    includelist = [];
end


%dataq = burst_peak(dataq);
[~, ~, ~, dataq] = proc_analysis(dataq, 1, 0, 3, 0);

%if there are manual corrections to burst2peak field, swap those in:
if ~isempty(intersect(construct, [9,16, 24, 14, 17, 15]))
     disp(['HP ' num2str(construct) ' data: swapping in peak position corrections.'])
    for i=1:numel(dataq)
        if isfield(dataq{i}, 'burst2peak2')
            dataq{i}.burst2peak_init = dataq{i}.burst2peak;
            dataq{i}.burst2peak = dataq{i}.burst2peak2;
            
        end
    end
end

end