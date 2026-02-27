%240906, AVT

% Calculates lifetimes of different constructs, subtracting out lifetimes of 
% a reference construct.
% Some relevant parameters are set in the HPlifetimes_errors code!

% 220804 note: ranges separate in physical space, but adjacent in HP7 space,
% appear continuous. Rebinning is done before conversion and thus bleeds
% from across the boundaries. Accept for now, maybe refine later.

% 240205: updated binning to reflect earlier changes in HPlifetimes_errors.
% Since that code now returns bin_centers, don't need to shift by half a
% bin in plots here.

% 250605: updated to include burst list with defined start and end points
% (timed_burst_list).

% Sample usage:
% [lifetime_ratio, lifetime_ratio_errors] = HPlifetimes_errors_differences(0, 1, 25, 1, data);
% [lifetime_ratio, lifetime_ratio_errors] = HPlifetimes_errors_differences(2, 1, 25, 1, data, [9 15 14 16 17]);
%[~, ~] = HPlifetimes_errors_differences(0, 1, 25, 0, data, [19 27 23 22 20 21]);
%[~, ~] = HPlifetimes_errors_differences(2, -1, 25, 0, data, [19 20 21 22 23], rzlist);
%[~, ~] = HPlifetimes_errors_differences(0, -1, 25, 0, data, [19 20 21 22 23], rz_IncludeBurstList_btstrp_3set);
%[~, ~] = HPlifetimes_errors_differences(0, 1, 25, 0, data, [9 24 15 14 16 17], uw_IncludeBurstList_btstrp_5set);

function [lifetime_diff, lifetime_diff_errors]= HPlifetimes_errors_differences(rebinhw, uw_only, proc_cutoff, ...
    baselineselect, data, construct, burst_list, timed_burst_list)

global HPcolours
global HPnames

plot_type = 'Mean';
%plot_type = 'Median';

binw=1;
hp7space =0;
%construct=7:11;
%construct = [9 15 14 16];
%construct = [7 8 9];

 % reference construct: calculating ratios of lifetimes with respect to it
refc = construct(1);
%refc = 9;

if hp7space
    binw =1;
    disp('Using bin width 1, for correct conversion to hp7 space.')
    %bins = 1:binw:76;
    %bins=0:binw:76; not needed -- get bin_centers out of HPlifetimes_errors.
%else
    %bins = 0:binw:(binw*floor(76/binw)); not needed -- get bin_centers out of HPlifetimes_errors.
end

if nargin <6
    construct=[19 20 21 22 23];
end

if nargin<7 || isempty(burst_list)
   burst_list =cell(1, max(construct)); % empty placeholder --> all bursts get included
end

if nargin<8 || isempty(burst_list)
   timed_burst_list =cell(1, max(construct)); % empty placeholder --> custom-range analysis skipped
end

if nargin <3       
    disp(['Loading data for constructs ' (num2str(construct))])
    for c=construct
        data{c} =  DatasetLoading_Baselineselect(c, baselineselect);
    end
else
    disp('Note: when dataset passed in as argument, baselineselect is not implemented.')
end

[mean_lifetime, sem, position_hists, bin_centers, avgforces] = HPlifetimes_errors(rebinhw, uw_only, 0,...
    proc_cutoff, plot_type,data, construct, burst_list, timed_burst_list);

% calculate differences of lifetimes

for c=construct
    
    if hp7space && sum(ismember([7:11 14 15], c))
        mean_lifetime{c}=Convert2HP7Space(mean_lifetime{c}, c, 0);
        sem{c}=Convert2HP7Space(sem{c}, c, 0);
    end
    
    lifetime_diff{c} = mean_lifetime{c}-mean_lifetime{refc};
    lifetime_diff_errors{c} = sqrt(sem{c}.^2 + sem{refc}.^2);
    
end


figure; hold on;


legendstr={};
    
for c=construct
    %default hairpin colours
    colors(c, :)= HPcolours{c};
    %%errorbar((bins(1:end-1)+binw/2),lifetime_ratio{c}, lifetime_ratio_errors{c},'Linestyle', 'none', 'Marker', 'o', 'Color', colors(c, :));
    %errorbar((bin_centers),lifetime_ratio{c}, lifetime_ratio_errors{c},'Linestyle', 'none', 'Marker', 'o', 'Color', colors(c, :));
    
    % instead of error bars, draw shaded regions
    if isequal(c, 16)
        %curtail CPD construct plotting where too few traces
%plotbins = find()
plot(bin_centers(1:44),lifetime_diff{c}(1:44), 'Color', colors(c, :));
        [px,py] = shadowerrorbar(bin_centers(1:44),lifetime_diff{c}(1:44),...
            lifetime_diff_errors{c}(1:44), lifetime_diff_errors{c}(1:44));
        fill(px, py, colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'None');

    else
        plot(bin_centers,lifetime_diff{c}, 'Color', colors(c, :));
        [px,py] = shadowerrorbar(bin_centers,lifetime_diff{c},...
            lifetime_diff_errors{c}, lifetime_diff_errors{c});
        fill(px, py, colors(c, :), 'FaceAlpha', 0.2, 'EdgeColor', 'None');
        %legendstr = [legendstr {['HP ' num2str(c)]}];
    end
    legendstr = [legendstr HPnames{c}];
end


if hp7space
    xlabel('Hairpin position, HP7 space (bp)'); 
else
 xlabel('Hairpin position (bp)'); 
end

 ylabel('Excess dwell time (s)'); 
 titlestr=[{[plot_type ' time at each position, wrt. ' HPnames{refc}]}, ...
     {['Rebinning ' num2str((1+rebinhw*2)) ' bp']}];
 if hp7space
    titlestr = [titlestr, {'HP7 space where applicable'}]; 
 end
 title(titlestr);
% title(['Ratio of mean time spent at each position, wrt constuct ' num2str(refc) newline 'HP7 space where applicable']);
xlims = xlim();ylims = (ylim);
% plot(xlims, [1 1], '--k');
plot([35 35], ylims, '--k');
 xlim([0 75]);
 f=get(gca, 'Children');
 legend(flipud(f(3:2:end)),legendstr); %legend('off')

end