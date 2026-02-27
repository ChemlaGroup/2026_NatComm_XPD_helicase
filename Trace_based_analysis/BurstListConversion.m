% 241107, AVT
% Standalone code to convert from burst_list, unequivocally defining
% bursts to include in an analysis, to indeces of the dataset in question.

% data: cell array corresponding to single dataset

% burst_list has the following in each row:
% date trace_number burst_number
% where all are doubles.

%Returns:
% -trace_num: vector of trace indeces within dataset for selected traces
% -burst_num: cell array with ith cell corresponding to ith trace in dataset;
%   relevant cells contain vector of selected burst indeces within trace.


function [trace_num, burst_num]=BurstListConversion(data, burst_list)

if isempty(burst_list)
    trace_num= (1:numel(data));
    burst_num = {}; % to include all bursts: burst number 0
    
else
    % select traces, bursts from subset list
    tlist = []; % for ease of indexing, make list in same format as burst_list from data structure
    trace_num = [];
    for i=1:numel(data)
        % separate entry for each trace-burst combination
        %tlist = [tlist; repmat([str2double(data{i}.date) data{i}.refnums(4)], [size(data{i}.burst2, 2), 1])];
        % single entry for each trace
        tlist = [tlist; [str2double(data{i}.date) data{i}.refnums(4)]];
    end
    
    burst_num = {};
    
    for i=1:size(burst_list, 1)
        tnum=find(tlist(:, 2)==burst_list(i, 2) & tlist(:, 1) ==burst_list(i, 1));
        if ~isempty(tnum)
            trace_num = [trace_num tnum]; % index of trace
            %burst_num{tnum}= [burst_num{tnum} burst_list(i, 3)] % index of burst within trace
        else
            error(['Cannot find trace with date ' num2str(burst_list(i, 1)) ' and trace ' num2str(burst_list(i, 2))]);
        end
    end
    
    trace_num=unique(trace_num);
    
    for i=trace_num
        inds= find(burst_list(:, 2)==data{i}.refnums(4) & abs(burst_list(:, 1)-str2double(data{i}.date))<10*eps);
        burst_num{i}(1, :) = burst_list(inds, 3); % index of burst within trace
        if size(burst_list(inds, :), 2) == 5
            burst_num{i}(2:3, :)= burst_list(inds, 4:5)';
            
        end
    end
    
end

end