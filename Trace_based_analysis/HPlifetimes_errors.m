% 220511, AVT

%Plots lifetimes with errorbars for all sequences in a set

% Returns a cell array (each element corresponds to one sequence) containing average time ('lifetime'/'residence time')
% spent at each position, calculated over the unwinding portions (optional) of all bursts, as well as the corresponding s.e.m.

% Inputs:

% -- rebinhw: allows rebinning of data (running window of half-width rebinhw) to minimize effect of binning artefacts
% (number of bins to use in rebinning, not basepairs)
% -- uw_only: which part of burst to use for calculations
%    - 0: whole burst
%    - 1: unwinding only
%    - -1: rezipping only
%    - 0.5: custom range (only if provided through timed_burst_list)
% -- makeplot: whether to display plot of lifetimes and errors
% -- proc_cutoff: only bursts reaching this position or beyond are included in analysis
% -- plot_type: 'Mean' or 'Median' of distribution returned at each hairpin
% position, with appropriate error
% -- baselineselect: whether to use quality cutoffs when loading dataset
% -- data: pre-loaded dataset (optional) -- required from 230918
% -- burst_list{c}{i} has the form:
%       date trace_number burst_number
% -- timed_burst_list{c}{i} has the form:
%       date trace_number burst_number start_time end_time

% If hist_scaled enabled: normalize the lifetimes and by average value for
% construct, recalculate s.e.m. Only affects plotting, not output.
% (Scaling: fudge factor at the moment -- not sure what causes the difference between datasets.)


%221031: changed bins to start at 0
%221230: updated title to reflect analysis parameters
%230212: updated error calculations: if error is zero (i.e. from single
%point), replace with NaN.
%231018: updated bins to be offset from whole hairpin positions by half a
%bin width; ensures histogram rounds to nearest value instead of taking
%floor.

%250526: updated to enable custom time-range selection.

%Sample usage:
%[~, ~, ~, ~] = HPlifetimes_errors(2, 1, 1, 25, 'Mean', data, [19 20 21 22 23]);
%[m_hist, sem, ~, ~] = HPlifetimes_errors(0, 1, 1, 25, 'Mean', data);
% [m_hist, sem, position_hists, bin_centers, avgforces]=HPlifetimes_errors(2, 1, 1, 25, 'Mean', data, [9 15 14 16 17 18]);
%[m_hist, sem, position_hists, bin_centers, avgforces]=HPlifetimes_errors(2, 1, 1, 25, 'Mean', data, [19 20 21], rzlist);
%[~, ~, ~, ~] = HPlifetimes_errors(0, 1, 1, 25, 'Mean', data, [9 24 15 14 16 17], uw_IncludeBurstList_btstrp_5set);

function [lifetime_hist, lifetime_errors, position_hists, bin_centers, avgforces] = HPlifetimes_errors(rebinhw, uw_only, makeplot, ...
    proc_cutoff, plot_type, data, construct, burst_list, timed_burst_list)

if nargin < 7
%construct = 7:11;
construct = [9 15 14 16];
%construct = [9 19 20];
%construct = [19 20 21];
%construct = [15 18];
%construct = [9];
%construct = [9 16 160171 90171];
%construct = [1]; %running on test dataset
end

if nargin<2
    uw_only=1;
    
end

if nargin<3
    makeplot=1;
end

if nargin<4
    proc_cutoff=25;
end

if nargin<5
    plot_type = 'Mean';
end

if nargin<6
    error('No data passed in!')
end

if nargin<8 || isempty(burst_list)
   burst_list =cell(1, max(construct)); % empty placeholder --> all bursts get included
end

if nargin<9 || isempty(burst_list)
   timed_burst_list =cell(1, max(construct)); % empty placeholder --> custom-range analysis skipped
end

% if nargin <7       
%     disp(['Loading data for constructs ' (num2str(construct))])
%     for c=construct
%         data{c} =  DatasetLoading_Baselineselect(c, baselineselect);
%     end
% else
%     disp('Note: when dataset passed in as argument, baselineselect value is ignored.')
% end

if uw_only == 1
   disp('Unwinding only!'); 
   sum_range= [20 50]; %[38 52]; %[18 50]; %[18 60]; % range for calculating sum of lifetimes
elseif uw_only == -1
    disp('Rezipping only! (All bursts, or subset?)'); sum_range=[20 50]%[18 42];%[18 34]%[43 53]; %[18 42]; %sum_range = [12 52]; %[20 50]; %rezipping
elseif uw_only == 0
    disp('Unwinding and rezipping!'); 
    sum_range= [38 52];
elseif uw_only == 0.5
    disp('Custom time range!'); sum_range=[18 42];
else 
    error('Please check uw_only value: should be 1, -1, or 0.')
end

binw=1;disp(['Bin width ' num2str(binw) ' bp']);
hp7space =0;
hist_scaled = 0; % only affects plotting; does not affect output
median_filtering=0;

if median_filtering
    filt_window = 0.3;
    disp('Using median filtering!')
else
    filt_window = [];
end

rebin_sum = 0; % if 1, rebinning returns sum over running window; if zero, returns average (directly comparable to non-rebinned). 
%binw = 1;


if hp7space
    binw =1;
    disp('Using bin width 1, for correct conversion to hp7 space.')
    %bins = 1:binw:76;
    bin_centers=0:binw:75;
else
    %bins=1:76;  % initial bins, before rebinning
    bin_centers = 0:binw:(binw*floor(75/binw));
end


%bins = (0-binw/2):binw:(binw*floor(75/binw)+binw/2); % calculating edges
lbinlim = 1+rebinhw; %number of bins
ubinlim = numel(bin_centers)-rebinhw-1; %number of bins
%rebin_inds = lbinlim:ubinlim;
%sum_range= [38 52]; %[18 50]; %[18 60]; % range for calculating sum of lifetimes
    
% Building histograms of binned position, then rebinning (running window of
% half-width rebinhw) to minimize effect of binning artefacts

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
   
    
    % Rebinning
    position_hists_rb{c}=position_hists{c}*NaN;  %edges padded with NaNs
    for b= lbinlim:ubinlim
        if rebin_sum
            % taking sum of rebinned set
            position_hists_rb{c}(:, b)=nansum(position_hists{c}(:, (b-rebinhw):(b+rebinhw)), 2); %disp('Rebinning with sum, not mean');
        else
            %taking mean of rebinning set
            position_hists_rb{c}(:, b)=nanmean(position_hists{c}(:, (b-rebinhw):(b+rebinhw)), 2);
        end
    end
    position_hists_rb{c}(isnan(position_hists{c}))=NaN;
    position_hists{c}=position_hists_rb{c};   
   
 
    if hp7space && sum(ismember([7:11 14 15], c))
        position_hists{c}=Convert2HP7Space(position_hists{c}, c, 0);
        disp('Converting to HP7 space!')
        %sem{c}=Convert2HP7Space(sem{c}, c, 0);
    end   
    
end


% Calculating lifetimes, s.e.m.

for c=construct
    
    survival{c}=sum(~isnan(position_hists{c}));
    
    if strcmp(plot_type, 'Mean')
        lifetime_hist{c} = nanmean(position_hists{c});    
        lifetime_errors{c}=nanstd(position_hists{c})./sqrt(survival{c});  % sem
    elseif strcmp(plot_type, 'Median')
        lifetime_hist{c} = median(position_hists{c}, 'omitnan');
        % error for median: calculation of notches for box plot
        
        quartiles = quantile(position_hists{c}, [0.25, 0.5, 0.75]); %
        iqr = quartiles(3, :)-quartiles(1, :);%interquartile range
        lifetime_errors{c} = 1.57*iqr./sqrt(survival{c});
    else
        error('Unknown plot type. Please input ''Mean'' or ''Median''. ');
    end    
    
    lifetime_errors{c}(lifetime_errors{c} < 1e-10)=NaN; % If only one burst, errors not defined
    
end
    
    titlestr =  [plot_type ' time spent at each position'];
   if median_filtering
       titlestr = [titlestr ' (with median filtering)'];
   end  
   titlestr2 = ['Rebinning window ' num2str((rebinhw*2+1))  ' bp'];
   if rebin_sum
      titlestr2 = [titlestr2 '; sum over window']; 
   else
      titlestr2 = [titlestr2 '; averaging over window']; 
   end
   titlestr3 = {titlestr; titlestr2};

if makeplot
    make_Mplot(construct, bin_centers, survival, lifetime_hist, lifetime_errors, Nbursts, hist_scaled);
    title(titlestr3);
    if hp7space
        xlabel('Hairpin position, HP7 space (bp)');
    else
        xlabel('Hairpin position (bp)');
    end
    
    [~, ~]= HPLifetimes_integrate(lifetime_hist, lifetime_errors, bin_centers, sum_range, construct);
    title(['Cumulative ' plot_type ' time spent over positions ' num2str(sum_range) ' bp']);
else 
    disp(titlestr3);
end


end

function make_Mplot(construct, bin_centers, survival, lifetime_hist, lifetime_errors, Nbursts, hist_scaled)

global HPcolours
global HPnames
legendstr={};

colors = distinguishable_colors(numel(construct));  %default colours

%survival probabilities
legendstr1={};
figure; 
subplot(3, 1, [1 2]); hold on;
for c=construct
    if ~isempty(HPcolours{c}) %sum(ismember([7:12 14:16, 17:20, 90171, 160171], c))
        colors(c, :)= HPcolours{c}; % standard hairpin colour scheme
        names{c} = HPnames{c};
    else
        disp('No colour preset for construct, using default.')
        colors(c, :)=[0 0 0];
        names{c} = ['Construct ' num2str(c)];
    end
    %plot((bin_centers(1:end-1)+binw/2),survival{c}, 'Color', colors(c, :), 'DisplayName', names{c})
    plot(bin_centers,survival{c}, 'Color', colors(c, :), 'DisplayName', names{c})
    %legendstr1 = [legendstr1 {HPnames{c}}];
end
xlabel('Hairpin position (bp)'); ylabel('Surviving bursts'); legend; %legend(legendstr1);
title('Number of bursts contributing at each hairpin position');

subplot(3, 1, 3); hold on;
for c=construct
 plot((bin_centers(2:end)),diff(survival{c}), 'Color', colors(c, :));
end
plot((bin_centers(1:end-2)), zeros(1, (numel(bin_centers)-2)), '--k') 
xlabel('Hairpin position (bp)'); ylabel('Change in burst #'); %legend(legendstr1);
    ylims = ylim;
    ylim([ylims(1) 5]);

%lifetimes and errors
figure; hold on;
for c=construct
    % shadowplot version:
    plot(bin_centers,lifetime_hist{c}, 'Color', colors(c, :), 'DisplayName', [names{c} ', ' num2str(Nbursts(c)) ' bursts']);
        [px,py] = shadowerrorbar(bin_centers,lifetime_hist{c}, lifetime_errors{c}, lifetime_errors{c});
        fill(px, py, colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
    
    %earlier errorbar version:
    %errorbar(bin_centers, lifetime_hist{c}, lifetime_errors{c},'Linestyle', 'none',...
    %    'Marker', 'o', 'Color', colors(c, :), 'DisplayName', [names{c} ', ' num2str(Nbursts(c)) ' bursts']);
end
xlabel('Hairpin position (bp)'); ylabel('Lifetime (s)');
% set(gca, 'XGrid', 'on')
% set(gca, 'XMinorGrid', 'on')
% set(gca, 'YGrid', 'on')
% set(gca, 'YMinorGrid', 'on')

% Scaling several of the above quantities to compare the constructs more
% effectively -- affects plotting only.
%(fudge factor at the moment -- not sure what causes the difference between datasets)

if hist_scaled
    for c=construct
        overall_mean_hist{c}= mean(lifetime_hist{c}(:), 'omitnan');
        lifetime_hist{c} = lifetime_hist{c}./overall_mean_hist{c};
        %lifetime_hist{c} = nanmean(lifetime_hist{c});
        lifetime_errors{c}=lifetime_errors{c}./overall_mean_hist{c};
    end
    for c=construct
        
  % shadowplot version:
    plot(bin_centers,lifetime_hist{c}, 'Color', colors(c, :), 'DisplayName', [names{c} ', ' num2str(Nbursts(c)) ' bursts; scaled']);
        [px,py] = shadowerrorbar(bin_centers,lifetime_hist{c}, lifetime_errors{c}, lifetime_errors{c});
        fill(px, py, colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
        
        %earlier errorbar variant:
       % errorbar(bin_centers,lifetime_hist{c}, lifetime_errors{c}, 'Linestyle', 'none', ...
        %    'Marker', 'o', 'Color', colors(c, :), 'DisplayName', [names{c} ', ' num2str(Nbursts(c)) ' bursts; scaled']);

    end

    ylabel('Lifetime (s)/ relative to mean');
end


end




