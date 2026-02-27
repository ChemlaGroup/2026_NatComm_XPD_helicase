% 231026, AVT

% Based on HPlifetimes_errors.

% Calculates the lifetime histogram (i.e. for all hairpin positions) for each set of traces from a given
% construct; then bootstraps over the histograms to calculate mean/median
% at each hairpin position

% Updated to include 
%a) subset of bursts to be analyzed fully
%b) separate subset of bursts to be analyzed over custom time range
%(following HPlifetimes_errors, HPbinpositions_medianfilt)

% Previous sample usage:
%%[lifetime_histBS, errorsBS, position_histsBS, binsBS,avgforcesBS]=
% HPlifetimes_bootstrap(2, 1, 1, 25, 'Mean', data, [19 20 21]);

% Sample usage:
% HPlifetimes_bootstrap(2, 1, 1, 25, 'Mean', data, [19 20 21]);
% HPlifetimes_bootstrap(0, 1, 1, 25, 'Mean', data, [9]);

function HPlifetimes_bootstrap(rebinhw, uw_only, makeplot, proc_cutoff, plot_type, data, construct, burst_list, timed_burst_list)

global HPnames


if nargin < 7
    %construct = 7:11;
    %construct = [9 15 14 16 17 18];
    construct = [19 20 21];
    %construct = [9 16 160171 90171];
    %construct = [1]; %running on test dataset
end

if nargin<2; uw_only=1; end

if nargin<3; makeplot=1;end

if nargin<4; proc_cutoff=25; end

if nargin<5;plot_type = 'Mean'; end

if nargin<6; error('No data passed in!');end

if nargin<8 || isempty(burst_list)
   burst_list =cell(1, max(construct)); % empty placeholder --> all bursts get included
end

if nargin<9 || isempty(timed_burst_list)
   timed_burst_list =cell(1, max(construct)); % empty placeholder --> custom-range analysis skipped
end

if strcmp(plot_type, 'Mean')
    bootfun = @mean_lifetime;
elseif strcmp(plot_type, 'Median')
    bootfun = @median_lifetime;
else
    error('Input statistic not implemented.')
end

if uw_only == 1
   disp('Unwinding only!'); direction = 'unwinding';
   sum_range= [18 40]; %[38 52]; %[18 50]; %[18 60]; % range for calculating sum of lifetimes
elseif   uw_only == -1
    disp('Rezipping only! (All bursts, or subset?)'); direction = 'rezipping';
    sum_range=[18 42]; %sum_range = [12 52]; %[20 50]; 
elseif uw_only == 0
    disp('Unwinding and rezipping!');  direction = 'unwinding and rezipping';
    sum_range= [38 52];
else 
    error('Please check uw_only value: should be 1, -1, or 0.')
end

binw=1;disp(['Bin width ' num2str(binw) ' bp']);

Niter = 30000;

hist_scaled = 0;
median_filtering=0;
if median_filtering
    filt_window = 0.3;
    disp('Using median filtering!')
else
    filt_window = [];
end


rebin_sum = 0; % if 1, rebinning returns sum over running window; if zero, returns average (directly comparable to non-rebinned). 
if rebin_sum
    disp('Rebinning with sum, not mean');
end
bin_centers = 0:binw:(binw*floor(76/binw));
BinningIndexOffset = 1-bin_centers(1);
lbinlim_ind = bin_centers(1)+rebinhw+BinningIndexOffset; % index of first rebinning bin (center)
ubinlim_ind = bin_centers(end)-rebinhw+BinningIndexOffset; % center of last rebinning bin (center)
%rebin_inds = lbinlim:ubinlim;
%sum_range = [18 50]; %[18 60]; % range for calculating sum of lifetimes

% Building histograms of binned position; rebinning before
% bootstrapping.

for c=construct
    
    if isempty(timed_burst_list{c}) % either all bursts or burst_list selection
        [position_hists{c}, Nbursts(c), avgforces{c}] = HPbinpositions_medianfilt(data{c}, bin_centers, uw_only, proc_cutoff, filt_window, burst_list{c});
    elseif ~isempty(burst_list{c}) && ~isempty(timed_burst_list{c})
        % burst_list and timed_burst_list selections: do both and combine
        [position_hists_1, Nbursts_1, avgforces_1] = HPbinpositions_medianfilt(data{c}, bin_centers, uw_only, proc_cutoff, filt_window, burst_list{c});
        [position_hists_2, Nbursts_2, avgforces_2] = HPbinpositions_medianfilt(data{c}, bin_centers, 0.5, proc_cutoff, filt_window, timed_burst_list{c});
        position_hists{c}= vertcat(position_hists_1, position_hists_2);
        Nbursts(c)=Nbursts_1 + Nbursts_2;
        avgforces{c}=horzcat(avgforces_1, avgforces_2);
    elseif ~isempty(timed_burst_list{c}) % only timed_burst_list selection
        [position_hists{c}, Nbursts(c), avgforces{c}] = HPbinpositions_medianfilt(data{c}, bin_centers, 0.5, proc_cutoff, filt_window, timed_burst_list{c});
    end
    
    
   
    if isempty(position_hists{c})
        error(['No suitable data for construct ' num2str(c) '.'])
    end

    
% Rebinning (running window of half-width rebinhw) to minimize effect of binning artefacts
     position_hists_rb{c}=position_hists{c}*NaN;  %edges padded with NaNs
    for b= lbinlim_ind:ubinlim_ind
        if rebin_sum
            % taking sum of rebinned set
            position_hists_rb{c}(:, b)=nansum(position_hists{c}(:, (b-rebinhw):(b+rebinhw)), 2); 
        else
            %taking mean of rebinning set
            position_hists_rb{c}(:, b)=nanmean(position_hists{c}(:, (b-rebinhw):(b+rebinhw)), 2);
        end
    end
    position_hists_rb{c}(isnan(position_hists{c}))=NaN;
    % replace position_hists with rebinned version
    position_hists{c}=position_hists_rb{c};
    
    %survival{c}=sum(~isnan(position_hists{c})); % do I need this?
end


% Bootstrap over position histograms to get mean lifetimes

for c=construct

    % replace NaN columns with -1 (deliberately unphysical) for
    % bootstrapping, then switch to NaNs again after:
    %nan_inds = find(isnan(sum(position_hists{c}, 'omitnan')));
    %position_hists{c}(nan_inds)=-1;
    position_hists{c}(:, 1:rebinhw) = -1;
    position_hists{c}(:,(end-rebinhw+1):end)=-1;
    
    tic
    %[ci, NKDtot] = bootci(Niter, bootfunNorm, bootstrap_input, max_pos_input, step, renorm_lowerbound, renorm_upperbound, 0, renormalize_NKDs);
    [ci{c}, lifetimes_statistic{c}]=bootci(Niter, bootfun, position_hists{c}, rebinhw);
    toc
    
    position_hists{c}(:, 1:rebinhw) = NaN;
    position_hists{c}(:,(end-rebinhw+1):end)=NaN;
    ci{c}(1:rebinhw) = NaN;ci{c}((end-rebinhw+1):end)=NaN;
    lifetimes_statistic{c}(:, 1:rebinhw) = NaN;lifetimes_statistic{c}(:, (end-rebinhw+1):end)=NaN;
   
    %lifetime_stat_std{c} = std(lifetimes_statistic{c});
end
% plot distribution of lifetimes statistic from bootstrapping



if makeplot
    %figure; hold on;
    for c = construct
         %errorbar(bin_centers, lifetimes_statistic{c}, ci{c},'Linestyle', 'none',...
        %'Marker', 'o', 'Color', colors(c, :)); %, 'DisplayName', [names{c} ', ' num2str(Nbursts(c)) ' bursts']);
    %legendstr = [legendstr {['HP ' num2str(c) ', ' num2str(Nbursts(c)) ' bursts']}];
        [h, handles]=Bootstrap_hppos_plots(lifetimes_statistic{c}, bin_centers, [1, 1, 1]);
        figure(handles(1));
        title({[plot_type ' bootstrapped lifetime, ' direction ', ' HPnames{c}]; ['Mean over ' num2str(rebinhw*2+1) ' bp bins, proc cutoff ' num2str(proc_cutoff)]});
           xlabel('Dwell time (s)');
        for j=2:3
           figure(handles(j));
           title({[plot_type ' bootstrapped lifetime, ' direction ', ' HPnames{c}]; ['Mean over ' num2str(rebinhw*2+1) ' bp bins, proc cutoff ' num2str(proc_cutoff)]});
           ylabel('Dwell time (s)');
        end
        
        % plot also original lifetime distribution overtop
    orig_stat{c}=bootfun(position_hists{c}, rebinhw);
    orig_stat{c}(1:rebinhw) = NaN;orig_stat{c}((end-rebinhw+1):end)=NaN;
    figure(handles(end)); % contour plot
    plot(bin_centers,orig_stat{c}, 'k', 'DisplayName', 'Original distribution')
    plot(bin_centers, ci{c}(1, :), '--r')
    plot(bin_centers, ci{c}(2, :), '--r')
    end

  
end



end

function mean_pos_lifetime = mean_lifetime(position_hist_mat, rebinhw)

mean_pos_lifetime = mean(position_hist_mat(:, (rebinhw+1):(end-rebinhw)), 'omitnan');
mean_pos_lifetime = [ones(1, rebinhw) mean_pos_lifetime ones(1, rebinhw)];

end

function median_pos_lifetime = median_lifetime(position_hist_mat, rebinhw)

median_pos_lifetime = median(position_hist_mat(:, (rebinhw+1):(end-rebinhw)), 'omitnan');
median_pos_lifetime = [ones(1, rebinhw) median_pos_lifetime ones(1, rebinhw)];

end

