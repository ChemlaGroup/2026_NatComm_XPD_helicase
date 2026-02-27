
% Note: this code does not autonomously classify into all categories the
% instances where XPD interacts with the modification. It did preliminary
% sorting which was refined by manual selection.

%% Get matrix of burst processivities for all traces
% need data already loaded

c= 14; % construct

Ntraces{c}=numel(data{c});

%downsample fleezers data -- test
% for i= 1:Ntraces{c}
%     data{c}{i}.bp= downsample(data{c}{i}.bp, 7);
%     data{c}{i}.time= downsample(data{c}{i}.time, 7);
%     %data{c}{i}.bp= decimate(data{c}{i}.bp, 7);
%     %data{c}{i}.time= decimate(data{c}{i}.time, 7);
% end
% data{c} = burst_peak(data{c});
% [~, ~, ~, data{c}] = proc_analysis(data{c}, 1, 0, 3, 0);


maxNb=1;

for i= 1:Ntraces{c}
    maxNb = max(maxNb, numel(data{c}{i}.proc));
end

burst_procs = NaN(maxNb, Ntraces{c});


for j= 1:Ntraces{c}
    i= numel(data{c}{j}.proc);
    burst_procs(1:i, j) = data{c}{j}.proc';
end

%% find bursts of different categories

%reaching at least 35 bp (rounded): total bursts

%ideally, would set to 34.5, but because of baseline issues,
%modification-associated plateau sometimes appears lower.
mod_le = 30.5; % lower estimate
mod = 34.5;

[burstInds_tot, traceInds_tot] = find(burst_procs>mod_le);


% clearing the modification

clear_threshold = 15; % how far from the modification the fork should be for XPD considered to be clear of the modification

[burstInds_clear, traceInds_clear] = find(burst_procs>=(mod+clear_threshold));


% estimated type I behaviour: unwinding does not progress beyond the modification. Put in a few
% basepairs extra as margin for fluctuation. Will have some overlap with type IV, dissociation.
margin = 3; 
[burstInds_block, traceInds_block] = find(burst_procs>mod_le & burst_procs<(mod+margin));

% estimated type II/III/IV behaviour: moving beyond modification, but not
% enough to clear; then either (II) backsliding, (III) rezipping, or (IV) dissociating.

%offset = 10;
[burstInds_adv, traceInds_adv] = find(burst_procs>=(mod+margin) & burst_procs<(mod+clear_threshold));

 % everything else:
%traceInds_class = [traceInds_clear, traceInds_block, traceInds_adv];
%burstInds_class = [burstInds_clear, burstInds_block, burstInds_adv];
%leftoverInds = setxor([traceInds_tot' burstInds_tot'], [traceInds_class' burstInds_class'], 'rows');
 
%% kernel density of burst processivities

%select subset
BurstInds = burstInds_tot; TraceInds= traceInds_tot;

figure; hold on;
step_err = 1*ones(1, numel(BurstInds));
[~, ~, kern] = kernel_density(burst_procs(sub2ind([maxNb, Ntraces{c}], BurstInds, TraceInds)), step_err);
title(['Kernel density of burst processivity (\geq 31 bp ) for ' HPnames{c}])

%% plotting all bursts reaching modification

BurstInds = burstInds_tot;
TraceInds= traceInds_tot;

figure; hold on;
k=numel(TraceInds);
plottedTraces = [];

for i=1:k
    
    t0 = data{c}{TraceInds(i)}.burst2(1,1);
    if isempty(intersect(TraceInds(i), plottedTraces))
        plot((data{c}{TraceInds(i)}.time-t0), data{c}{TraceInds(i)}.bp, 'k', ...
            'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4))])
    end
    
    t1 = data{c}{TraceInds(i)}.burst2(1, BurstInds(i));
    t2 = data{c}{TraceInds(i)}.burst2(2, BurstInds(i));
    plot_inds = find(data{c}{TraceInds(i)}.time >t1 & data{c}{TraceInds(i)}.time<t2);
    plot((data{c}{TraceInds(i)}.time(plot_inds)-t0), data{c}{TraceInds(i)}.bp(plot_inds), 'b',...
        'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4)) ' ' num2str(BurstInds(i))])
    plottedTraces = [plottedTraces TraceInds(i)];
        
end

plot([0 100], [35 35], 'k:', 'DisplayName', 'Modification')


xlabel('Time from start of first burst (s)')
ylabel('Basepairs unwound (bp)')
ylim([-20 100]);
title('Bursts reaching modification')

xlim([-5 100])
zoom reset


%% Figures to confirm/identify burst classification

% bursts clearing modification
BurstInds = burstInds_clear;
TraceInds= traceInds_clear;

figure; hold on;
k=numel(TraceInds);
plottedTraces = [];

for i=1:k
    
    t0 = data{c}{TraceInds(i)}.burst2(1,1);
    if isempty(intersect(TraceInds(i), plottedTraces))
        plot((data{c}{TraceInds(i)}.time-t0), data{c}{TraceInds(i)}.bp, 'k', ...
            'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4))])
    end
    
    t1 = data{c}{TraceInds(i)}.burst2(1, BurstInds(i));
    t2 = data{c}{TraceInds(i)}.burst2(2, BurstInds(i));
    plot_inds = find(data{c}{TraceInds(i)}.time >t1 & data{c}{TraceInds(i)}.time<t2);
    plot((data{c}{TraceInds(i)}.time(plot_inds)-t0), data{c}{TraceInds(i)}.bp(plot_inds), 'g',...
        'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4)) ' ' num2str(BurstInds(i))])
    plottedTraces = [plottedTraces TraceInds(i)];
        
end

plot([0 100], [35 35], 'k:', 'DisplayName', 'Modification')
plot([0 100], [(mod+clear_threshold) (mod+clear_threshold)], 'k:', 'DisplayName', 'Cleared')

xlabel('Time from start of first burst (s)')
ylabel('Basepairs unwound (bp)')
ylim([-20 100]);
title('Bursts clearing modification')

xlim([-5 100])
zoom reset


%%

% bursts blocked by modification
BurstInds = burstInds_block;
TraceInds= traceInds_block;

figure; hold on;
k=numel(TraceInds);
plottedTraces = [];

for i=1:k
    
    t0 = data{c}{TraceInds(i)}.burst2(1,1);
    if isempty(intersect(TraceInds(i), plottedTraces))
        plot((data{c}{TraceInds(i)}.time-t0), data{c}{TraceInds(i)}.bp, 'k', ...
            'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4))])
    end
    
    t1 = data{c}{TraceInds(i)}.burst2(1, BurstInds(i));
    t2 = data{c}{TraceInds(i)}.burst2(2, BurstInds(i));
    plot_inds = find(data{c}{TraceInds(i)}.time >t1 & data{c}{TraceInds(i)}.time<t2);
    plot((data{c}{TraceInds(i)}.time(plot_inds)-t0), data{c}{TraceInds(i)}.bp(plot_inds), 'r',...
        'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4)) ' ' num2str(BurstInds(i))])
    plottedTraces = [plottedTraces TraceInds(i)];
        
end

plot([0 100], [35 35], 'k:', 'DisplayName', 'Modification')

xlabel('Time from start of first burst (s)')
ylabel('Basepairs unwound (bp)')
ylim([-20 100]);
title('Bursts with unwinding blocked by modification (and/or excursions)')
ylim([-10 45])
xlim([-5 100])
zoom reset

%%

% bursts advancing only partially beyond modification
BurstInds = burstInds_adv;
TraceInds= traceInds_adv;

figure; hold on;
k=numel(TraceInds);
plottedTraces = [];

for i=1:k
    
    t0 = data{c}{TraceInds(i)}.burst2(1,1);
    if isempty(intersect(TraceInds(i), plottedTraces))
        plot((data{c}{TraceInds(i)}.time-t0), data{c}{TraceInds(i)}.bp, 'k', ...
            'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4))])
    end
    
    t1 = data{c}{TraceInds(i)}.burst2(1, BurstInds(i));
    t2 = data{c}{TraceInds(i)}.burst2(2, BurstInds(i));
    plot_inds = find(data{c}{TraceInds(i)}.time >t1 & data{c}{TraceInds(i)}.time<t2);
    plot((data{c}{TraceInds(i)}.time(plot_inds)-t0), data{c}{TraceInds(i)}.bp(plot_inds), 'Color', [1 0.4 0.15],...
        'DisplayName', [data{c}{TraceInds(i)}.date ' ' num2str(data{c}{TraceInds(i)}.refnums(4)) ' ' num2str(BurstInds(i))])
    plottedTraces = [plottedTraces TraceInds(i)];
        
end

plot([0 140], [35 35], 'k:', 'DisplayName', 'Modification')
plot([0 140], [(mod+clear_threshold) (mod+clear_threshold)], 'k:', 'DisplayName', 'Cleared')
xlabel('Time from start of first burst (s)')
ylabel('Basepairs unwound (bp)')
ylim([-20 100]);
title('Bursts partially advancing beyond modification')

ylim([-10 55])
xlim([-5 100])
zoom reset

%% Printing index list

datelist = [];
trace_nums = [];

for i=1:numel(TraceInds)
    datelist = [datelist str2num(data{c}{TraceInds(i)}.date)];
trace_nums = [trace_nums data{c}{TraceInds(i)}.refnums(4)];
end

[datelist' TraceInds trace_nums' BurstInds]

%% plotting burst subsets already sorted

% blist = [TraceInds trace_nums' BurstInds] % formatting only; load from
% spreadsheet

c=15;
mod = 34.5;
clear_threshold = 15;

figure; hold on;

for i=1:size(blist, 1)
    
    t0= data{c}{blist(i, 1)}.burst2(1, blist(i, 3));
    t_end = data{c}{blist(i, 1)}.burst2(2, blist(i, 3));
plot_inds = find(data{c}{blist(i, 1)}.time >=t0 & data{c}{blist(i, 1)}.time < t_end);
    plot((data{c}{blist(i, 1)}.time(plot_inds)-t0), data{c}{blist(i, 1)}.bp(plot_inds), 'Color', HPcolours{c},... %[1 0.4 0.15],...
        'DisplayName', [data{c}{blist(i, 1)}.date ' ' num2str(blist(i, 2)) ' ' num2str(blist(i, 3))])
end

plot([0 30], [35 35], 'k:', 'DisplayName', 'Modification')
%plot([0 30], [(mod+clear_threshold) (mod+clear_threshold)], 'k:', 'DisplayName', 'Cleared')
xlabel('Time from start of burst (s)')
ylabel('Basepairs unwound (bp)')
ylim([-5 (mod+clear_threshold+10)]);
xlim([-2 24]);



