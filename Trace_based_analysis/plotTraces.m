
% plot activity portion of all traces of a dataset, with t=0 being the start of activity (first
% point of first burst) for all traces

%Sample usage:
% [data14, ~] = DatasetLoading_Baselineselect(14, 1);
% plotTraces(data14, 25);

function plotTraces(data, proc_cutoff, median_filter)

if nargin<2
    proc_cutoff = 0; %default
    %proc_cutoff = 32;
    median_filter =0;
end

if nargin<3
   median_filter = 0; 
end

figure; hold on;
Ntraces = numel(data);
Ncolors = distinguishable_colors(Ntraces);
legendstr = {};
Toffset = 5; % (s) how much time before activity to plot

% More plotting parameters

%median_filter = 0;
medfilt_window = 0.3; %s
plot_baseline3 =1; % plot baseline after end of activity

for i=1:Ntraces
    if plot_baseline3
        activity_inds = find(data{i}.time >= (data{i}.burst2(1, 1)-Toffset));
    else
        activity_inds = find(data{i}.time >= (data{i}.burst2(1, 1)-Toffset) & data{i}.time <= data{i}.burst2(2, end));
    end
    
    if median_filter
        medfilt_npoints = round(data{i}.traprate*medfilt_window);
        bp = movmedian(data{i}.bp, medfilt_npoints);
    else
        bp = data{i}.bp;
    end
    
    % processivity cutoff: should not be affected by what happens
    % before/after activity.
    inds = find(data{i}.time >= data{i}.burst2(1, 1) & data{i}.time <= data{i}.burst2(2, end));
    if max(bp(inds)) >= proc_cutoff
        t0 =  data{i}.time(activity_inds(1))+Toffset;
       
        plot(data{i}.time(activity_inds)-t0, bp(activity_inds), 'Color', Ncolors(i, :));
        legendstr = [legendstr [data{i}.date '\_' num2str(data{i}.refnums(4))]];
    end
end

xlabel('Time from start of activity (s)')
ylabel('Bp unwound')
legend(legendstr);
legend('off')
titlestr=['Traces of HP ' data{1}.construct(7:end) '; proc cutoff: ' num2str(proc_cutoff) ' bp'];

if median_filter
    titlestr = [titlestr '; median filter ' num2str(medfilt_window) ' s'];    
end

title(titlestr);

xlims = xlim;
%plot(xlims, [35 35], '--r');
%plot(xlims, [0 0], '--k');

plot([0 100], [35 35], '--k');
plot([0 100], [0 0], '--k');
ylim([-5 100]);
set(gca, 'YGrid', 'on');

end