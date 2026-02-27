%221102, AVT
% Slightly modified version of HPbinpositions. This one allows data to be
% median-filtered for the analysis; I want to see how much of a difference
% this makes.

% Given a dataset, bins each trace by position.
% Each row of resulting position_hists matrix corresponds to a single
% burst: total time (in seconds) spent at given position.
% Below and above actual unwinding for a given burst, row is filled with NaNs.
%FiltWindow: in s

% 230108: fixed bug: apply median filtering before processivity cutoff!
%230107: return vector of forces (average force of trace) corresponding to
%bursts.

% 221031: fixed indexing to work properly for binwidth other than 1.
% Also fixed bug that set top and bottom bin to NaN always, rather than when empty.

% 230209: burst_inds selection must come with LE sign, not L sign, to work
% consistently.

%231028: for more accurate binning, should have hairpin position as centre
%of bin (i.e. round up/down to nearest hp position instead of taking floor).
% Updated to take in explicitly bin centers, from which bin edges are
% calculated.

% 231121: added rezipping-only option (uw_only = -1).

% 240102: adding burst_list: set of bursts to include (presumably manually selected).
% (Section of code dealing with this not well-written but works.)
% burst_list is a 3-column matrix: 
% -- 1st column is date,
% -- 2nd column is the number of the trace (refnums(4) in the data
% structure)
% -- 3rd column is number of burst within the trace

%240908: fixed bug in burst_list implementation; taking date into account as
%well as trace number in last selection step to remove ambiguity.

% Sample usage:
% [position_hists, Nbursts, avgforces]=HPbinpositions_medianfilt(data7, 0:1:75 , 1, 25, 0.3);

function [position_hists, Nbursts, avgforces] = HPbinpositions_medianfilt(data, bin_centers, uw_only, proc_cutoff, FiltWindow, burst_list)

global HPnames
 
if nargin<2
    bin_centers = 1:1:75;
end

if nargin<3
    uw_only = 1;
end

if nargin<4
    proc_cutoff = 25;
end

if nargin<5
    FiltWindow = []; %s
end

if nargin<6
    burst_list = [];
end

if isequal(uw_only, 0)
    disp('Note: processivity cutoffs not enabled when rezipping included')
end


%disp('Note: if integer bins entered, histrogram takes floor values instead of rounding!')
binw = bin_centers(2)-bin_centers(1);
bin_edges= [bin_centers-binw/2 bin_centers(end)+binw/2];

% Calculate processivities -- changing to calculating them after data
% loading. Otherwise manual corrections will be overwritten.
% Replace 0s by 1s to output plots.
%[~, ~, ~, data]= proc_analysis(data, 1, 0, 3, 0, [], 0);% changed maxproc-dependent field to 0 % updated format 210415
%data = burst_peak(data);
position_hists = [];
Nbursts=0;
avgforces = [];
makeplot = 0;

% Selecting subset of bursts based on manually chosen list
    
if isempty(burst_list)
    % if no subset selection: use all data
    trace_num = (1:numel(data));
    burst_num = {}; % to include all bursts: burst number 0
    
else
    [trace_num, burst_num]=BurstListConversion(data, burst_list);
end

%per trace
%debugging purposes:
if makeplot
figure; hold on;
end

for i=trace_num
    
    if ~isempty(FiltWindow)
        medfilt_npoints = round(data{i}.traprate*FiltWindow);
    else
        medfilt_npoints = 1;
    end
    

    if isempty(burst_num)
        Nburst = 1:size(data{i}.burst2, 2);
    else
        Nburst = burst_num{i}(1, :);
    end
    
    %per burst
    for j= 1:numel(Nburst)
        % apply median filtering before processivity cutoff, but
        % after identifying region of interest
        if isequal(uw_only, 0)    % Both unwinding and rezipping of each burst
            burst_inds = find(data{i}.time >= data{i}.burst2(1, Nburst(j)) & data{i}.time <= data{i}.burst2(2, Nburst(j)));
            bp = movmedian(data{i}.bp(burst_inds), medfilt_npoints);
            %data_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time < data{i}.burst2(2, j));
            burst_hist = histcounts(bp, bin_edges);
            Nbursts=Nbursts+1;
            avgforces = [avgforces data{i}.avgforce];
        elseif isequal(uw_only,1)  % unwinding parts of bursts only, and only those bursts passing processivity cutoff
            burst_inds = find(data{i}.time >= data{i}.burst2(1, Nburst(j)) & data{i}.time <= data{i}.burst2peak(Nburst(j)));
            bp = movmedian(data{i}.bp(burst_inds), medfilt_npoints);
            if max(bp) >=proc_cutoff
                %data_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time < data{i}.burst2peak(j));
                burst_hist = histcounts(bp, bin_edges);
                Nbursts=Nbursts+1;
                avgforces = [avgforces data{i}.avgforce];
            else %ignore low-processivity bursts
                continue
            end
        elseif isequal(uw_only, -1) %rezipping
            burst_inds = find(data{i}.time <= data{i}.burst2(2, Nburst(j)) & data{i}.time >= data{i}.burst2peak(Nburst(j)));
            bp = movmedian(data{i}.bp(burst_inds), medfilt_npoints);
            if max(bp) >=proc_cutoff
                %data_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time < data{i}.burst2peak(j));
                burst_hist = histcounts(bp, bin_edges);
                Nbursts=Nbursts+1;
                avgforces = [avgforces data{i}.avgforce];
            else %ignore low-processivity bursts
                continue
            end
        elseif isequal(uw_only, 0.5) %manually-defined segment of burst
            t1= burst_num{i}(2, j);
            t2= burst_num{i}(3, j);
            burst_inds = find(data{i}.time >= t1 & data{i}.time <= t2);
            bp = movmedian(data{i}.bp(burst_inds), medfilt_npoints);
            if max(bp) >=proc_cutoff
                %data_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time < data{i}.burst2peak(j));
                burst_hist = histcounts(bp, bin_edges);
                Nbursts=Nbursts+1;
                avgforces = [avgforces data{i}.avgforce];
            else %ignore low-processivity bursts
                continue
            end
        end
        
         % for debugging purposes; if uncommenting, also uncomment figure
         % command above.
         if makeplot
             t= (0:(numel(bp)-1))/data{i}.traprate;
             plot(t, bp, 'DisplayName', [data{i}.date '_' num2str(data{i}.refnums(4)) ', burst ' num2str(Nburst(j))]);
             
             title(['Burst portions histogrammed for construct: ' data{i}.construct]); ylabel('Basepairs unwound (bp)'); xlabel('Time(s)');
             ylim([-5 95]);
         end
        
        % don't want zeros in bins above actual reach of burst to mess up the stats
        [~, max_bin_ind]=find(bin_edges> max(bp), 1);
        if max_bin_ind < numel(bin_centers)
            burst_hist(max_bin_ind:end)=  NaN;
        end
        % ditto for bursts not starting at baseline
        [~, min_bin_inds]=find(bin_edges< min(bp));
        if ~isempty(min_bin_inds) && min_bin_inds(end) >1
            burst_hist(1:(min_bin_inds(end)-1)) = NaN;
        end
        
        burst_hist = burst_hist./data{i}.traprate; % convert from timepoints to time (s)
        try
            position_hists = vertcat(position_hists, burst_hist);
        catch
            disp([i, j])
        end
        
       
        
    end
    
end
end