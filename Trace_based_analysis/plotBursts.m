% AVT


% Plots bursts of a dataset passing a given processivity cutoff.
% One subplot is hairpin position as a function of time, each burst starting
% at 0 s. The other subplot is a histogram of the bursts by position.  

% Modified significantly on 230108 to be analogous to plotTraces, allowing median filtering.

% 240825: fixed histogramming, with bin centers being whole basepairs,
% treating bin centers and bin edges explicitly.

%250522: added downsampling for fleezers data
%250621: disabled. Better to downsample data directly during processing.

% Sample usage
% plot_bursts(data15, 25, 1, 0);

function plotBursts(data, proc_cutoff, uw_only, median_filter)

if nargin<2
    proc_cutoff = 25;
end

if nargin<3
    uw_only = 0;
end

if nargin <4
    median_filter =0;
elseif median_filter
    medfilt_window = 0.05; %0.3; %s
end

figure; 
ax1=subplot(1, 5, [1 4]); hold on;
set(gca, 'YGrid', 'on');
ax2=subplot(1, 5, 5); hold on;
linkaxes([ax1 ax2], 'y');
N = distinguishable_colors(110);
count = 0;
%burst_list = [];

bin_centers=0:90;
binw = bin_centers(2)-bin_centers(1);
bin_edges= [bin_centers-binw/2 bin_centers(end)+binw/2];

legendstr={};

for i = 1:size(data, 2)
    
    for j=1:size(data{i}.burst2, 2)
        
        if abs(uw_only)<10*eps
            burst_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time <= data{i}.burst2(2, j));
        elseif abs(uw_only-1)<10*eps  % unwinding only
            burst_inds = find(data{i}.time >= data{i}.burst2(1, j) & data{i}.time <= data{i}.burst2peak(j));
        elseif abs(uw_only+1)<10*eps  % rezipping only
            burst_inds = find(data{i}.time >= data{i}.burst2peak(j) & data{i}.time <= data{i}.burst2(2, j));
        else
            error('Unwinding and/or rezipping?')
        end
        
                
%         if ~data{i}.oldtrap
%             %disp('Note: decimation of fleezers data currently has a shift issue.')
%             dfactor = ceil(data{i}.traprate/100);
%             bp = FilterAndDecimate(data{i}.bp(burst_inds), dfactor);
%             time = FilterAndDecimate(data{i}.time(burst_inds), dfactor);
%             %bp = decimate(data{i}.bp(burst_inds), dfactor);
%             %time = decimate(data{i}.time(burst_inds), dfactor);
%         else
            bp=data{i}.bp(burst_inds);
            time=data{i}.time(burst_inds);
             dfactor =1;
%         end
        
        
         %median filtering only includes burst itself, nothing from
         %adjacent timewindows.
        if median_filter
            medfilt_npoints = round(data{i}.traprate*medfilt_window/dfactor);
            bp = movmedian(bp, medfilt_npoints);
        end

        if max(bp) >=proc_cutoff
            count = count+1;
            %burst_list = vertcat(burst_list, [i j max(bp)]);
            t0 =  time(1);
            l1=plot(ax1, time-t0, bp, 'Color', N(count, :));
            legendstr = [legendstr {[ num2str(data{i}.date) '_' num2str(data{i}.refnums(4)) '_' num2str(j)]}];
            %h=histcounts(bp, bins, 'Normalization', 'probability');
            h=histcounts(bp, bin_edges)/data{i}.traprate;
            l2=plot(ax2, h, bin_centers, 'Color', N(count, :));
            % Link corresponding lines on the two subplots, so as to be
            % able to take out given burst from both simultaneously:
            hlink = linkprop([l1, l2],'Visible');
            setappdata(l1, 'Visibility', hlink);
            %setappdata(l2, 'Visibility',  hlink);
            
        end
    end
end

subplot(ax2);
ylim('auto');xlim('auto');
%xlabel('Fraction')
xlabel('Time (s)')
legend(legendstr);
legend('off');
xlims = xlim(gca);
plot(xlims, [35 35], '--k');
plot(xlims, [50 50], '--r');
set(gca, 'YGrid', 'on');

subplot(ax1);
ylim([bin_edges(1) bin_edges(end)]);
titlestr = [data{1}.construct(5:end) ', processivity > ' num2str(proc_cutoff) ' bp'];
if median_filter
    titlestr = [titlestr '; median filter ' num2str(medfilt_window) ' s'];    
end
title(titlestr);
ylabel('Hairpin position (bp)')
xlabel('Time (s)')
legend(legendstr);
legend('off');

% These are useful to me for looking at modification-associated pauses; can
% be omitted otherwise!
plot([0 30], [35 35], '--k');
plot([0 30], [50 50], '--r');

end