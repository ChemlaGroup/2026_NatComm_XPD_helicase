% From Barbara
% Analyze dip and burst processivity


% Modified 180828 to save burst processivities as a field in each element
% of 'data'. Need to save 'data' after running to retain processivities per
% element.

% 190328: Modifying histogram to plot bar in the middle of each bin, not the lowest index of each bin.
% 191118: Removing artificial 90bp boundary on values in m.
% 191118: Adding maxproc field; maximum processivity per molecule.
% Note: shift can currently only be enabled for some of the plots. Double-check before using. 
% 201106: Allowing force to be passed in.
% 230914: including proc2 calculations. These normally are done during
% burst selection, but need to be re-done for re-processed data.

% Takes in: 
% -- a data structure; 
% -- width of force bin considered; 
% -- option to draw plot of processivity of all forces combined;
% -- bin size for processivity plots;
% -- option to plot counts separately by force;
% -- option to plot maximum processivity per trace (all forces and separately by force)

% Outputs:
% proc: 
% i-th row: statistics over i-th force
%(i, 1); force; 
%(i, 2): number of processive bursts; 
%(i,3): total number of bursts;
%(i, 4): number of processive molecules;
% (i, 5): total number of molecules

% m: vector of maximum bp unwound per burst at given force, i.e. m{1} is vector at lowest force, etc
% force: array of forces
% data: data structure modified with fields calculated in code (data.proc, data.maxproc, data.tproc)

                                            
% Example usage:
% [proc, m, force, data] = proc_analysis(data, 1, 1, 3,0,  [], 1)

% Formerly
%[proc, m, max_mol_proc, force, data, fig_handles] = proc_analysis(data, forcebinw, drawplot_tot, binsize, drawcounts_force_sep, max_per_molecule)

function [proc, m, force, data] = proc_analysis(data, forcebinw, drawplot_tot, binsize, drawcounts_force_sep, force, max_per_molecule)

if nargin <3
    drawplot = 0;
    binsize = 0;
    drawcounts_force_sep = 0;
    force = [];
    max_per_molecule = 0;
end

if nargin <4
    binsize = 5;
    drawcounts_force_sep = 0;
     force = [];
    max_per_molecule = 0;
end

if nargin <5
    drawcounts_force_sep = 0;
     force = [];
    max_per_molecule = 0;
end

if nargin <6
     force = [];
    max_per_molecule = 0;
end

if nargin <7
    max_per_molecule = 0;
end

%shift = 5; % to predict where distribution of processivities will be for new sequence. This only affects the final plot.
shift = 0;
normhist = 1; % Normalize histogram if 1; plot counts if 0.
fig_handles = [];

%forcebinw=3;    % 'resolution' in force
%force= [10] 
%force= [10:forcebinw:14]; % [9:12]; %[0,2,5,10,50]; 

forces = [];
for i=1:size(data, 2)
    forces = [forces data{i}.avgforce];
end


if isempty(force)
    min_f = round(min(forces));
    max_f = round(max(forces));
    force= [min_f:forcebinw:max_f];
end

%force= [9:forcebinw:12];
%force= [10:forcebinw:12];

N = distinguishable_colors(length(force));

proc = [force', zeros(length(force),4)];    % i-th row: statistics over i-th force
                                            %(i, 1); force; (i, 2): number of processive bursts; 
                                            %(i,3): total number of bursts;
                                            %(i, 4): number of processive
                                            %molecules; (i, 5): total
                                            %number of molecules

m = cell(1,length(force));  % vector of maximum bp unwound per burst at given force, i.e. m{1} is vector at lowest force, etc.
mtracemax = cell(1,length(force));
mbp = []; err = [];
tavg = m;   % vector per force; time between start of first burst and crossing of processivity 'boundary' per molecule
mt = [];
l = [];
for i = 1:length(force)
    for k = 1:length(data)
         %if data{k}.conc{2,1}==conc(i)
         
         % calculate proc2, as it is not automatically calculated for
         % re-processed data:
         data{k}.proc2 = [];
         for b = 1:size(data{k}.burst2,2)
             ind = find(data{k}.time >= data{k}.burst2(1,b) & data{k}.time <= data{k}.burst2(2,b));
             
             if ~isempty(nonzeros((data{k}.bp(ind)>25)))
                 data{k}.proc2 = [data{k}.proc2, 1];
             else
                 data{k}.proc2 = [data{k}.proc2, 0];
             end
         end
         
        if abs(data{k}.avgforce - force(i)) <= forcebinw/2;          %round(data{k}.avgforce)==force(i) % changed to allow bins <1 in total width
            if ~isempty(data{k}.burst2)
                p = []; tp = [];   
                for j = 1:size(data{k}.burst2,2)   %iterating over each burst, (t_start, t_end)
                    
                    ind = find(data{k}.time > data{k}.burst2(1,j) & data{k}.time < data{k}.burst2(2,j));
                    bp = data{k}.bp(ind);
                    time = data{k}.time(ind); time = time - time(1); % time set to 0 at start of burst
                    if j ==1
                        t1 = data{k}.time(ind(1));
                    end
                    data{k}.proc(j)=max(data{k}.bp(ind));  % save processivity per burst
                    m{i} = [m{i} max(data{k}.bp(ind))];
                    if numel(data{k}.proc2)<j
                        disp([k j])
                    end
                    p = [p data{k}.proc2(j)];   % vector containing binary of processive bursts, per trace
                end
                
                data{k}.maxproc = max(data{k}.proc); % Maximum processivity per molecule
                %disp(data{k}.maxproc)
                mtracemax{i} = [mtracemax{i} data{k}.maxproc];
                
                proc(i,2) = proc(i,2)+sum(p); %sum(data{k}.proc2);
                proc(i,3) = proc(i,3)+length(p); %length(data{k}.proc2);
                
                
                if sum(p) >0 %sum(data{k}.proc2) >0
                    proc(i,4) = proc(i,4)+1;
                    p = find(data{k}.proc2>0);
                    ind = find(data{k}.time > data{k}.burst2(1,p(1)) & data{k}.time < data{k}.burst2(2,p(1)));
                    t = data{k}.time(ind); bp = data{k}.bp(ind);
                    t = t(bp>=25);
                    data{k}.tproc = t(1) - data{k}.burst2(1,1); % time between start of first burst and first crossing of processivity 'boundary' 
                     tavg{i} = [tavg{i} data{k}.tproc]; 
                else
                    
                    data{k}.tproc = inf;
                end
                proc(i,5) = proc(i,5)+1;
                
                
            end
        end
        %m
         %end
         
    end
    %m{i}(m{i}>90)=90; 
    l = [l length(m{i})];
    mbp(i) = mean(m{i}); err(i) = std(m{i})/sqrt(length(m{i}));
    mt(i) = mean(tavg{i}); terr(i) = std(tavg{i})/sqrt(length(tavg{i}));
end

if drawplot_tot
    
    
 figure; 
hold on;
plotSpread(m,'xValues',force,'spreadWidth',1.2,'distributionMarkers','o'); % length(conc)
errorbar(force,mbp,err,'o-k','LineWidth',2);
xlabel('Force (pN)'); ylabel('Processivity per burst (bp)');

%xlim([0 6]);

% figure;
% hold on;
% errorbar(force,mt,terr,'ko');
% title('Mean time to cross processivity boundary')
% xlabel('Force (pN)'); ylabel('Time (s)')
    
    h1 = figure; 
    fig_handles = [fig_handles h1]; 
    hold on;
    binned = cell(1,length(force));
    plotvec = zeros(length(0:binsize:100),1);   % total bursts, all forces combined
    lvec = [];

    for i=1:length(m)
        binned{i} = histc(m{i}+shift, (0:binsize:100));
        %plotvec = [plotvec binned{i}'];
        if ~isempty(binned{i})
            plotvec = plotvec + binned{i}';
        end
    end

    if normhist
        plotvec = plotvec/sum(plotvec);
    end
    
    plot((binsize/2:binsize:100+binsize/2), plotvec)

    xlabel('Max hairpin position reached (bp)','FontSize',16,'FontWeight','bold');
    
    if normhist 
        ylabel('Frequency','FontSize',16,'FontWeight','bold');
    else
        ylabel('Counts','FontSize',16,'FontWeight','bold');
    end
    xlim([0 100]);
    
    title(['Maximum hp position reached during burst, all forces; ' num2str(binsize) ' bp bins'], 'FontSize',16,'FontWeight','bold');

    
    
    
end

if drawcounts_force_sep
    h2 = figure; 
    fig_handles = [fig_handles h2];
    hold on;
    binned = cell(1,length(force));
    plotvec = zeros(length(0:binsize:100),length(force));
    
    
    for i=1:length(force)
        binned{i} = histc(m{i}+shift, (0:binsize:100));
        if isempty(binned{i})
            binned{i} = zeros(1,size(binned{i}, 1));
        end
        
        legendstr{i} = [num2str((force(i)-forcebinw/2), '% .1f') ' to ' num2str(force(i)+forcebinw/2, '% .1f') 'pN (' num2str(length(m{i})) ' bursts; ' num2str(proc(i, 5))  ' molecules)'];
        plot((binsize/2:binsize:100+binsize/2), binned{i}, 'Color', N(i, :));
        
    end
    
    xlim([0 100]);
    legend(legendstr)
    xlabel('Max hairpin position reached (bp)','FontSize',16,'FontWeight','bold');
    ylabel('Counts','FontSize',16,'FontWeight','bold');
    title(['Maximum hp position reached during burst at different forces, ' num2str(binsize) ' bp bins'], 'FontSize',16,'FontWeight','bold');
    
end

if max_per_molecule
    max_mol_proc = [];
    for k = 1:length(data)
        max_mol_proc = [max_mol_proc data{k}.maxproc];
    end
    
%     bin_mol_proc = histc(max_mol_proc, (0:binsize:100));
%     h3 = figure;
%     fig_handles = [fig_handles h3];
%     hold on;
%     plot((binsize/2:binsize:100+binsize/2), bin_mol_proc)
%     xlabel('Max hairpin position reached (bp)','FontSize',16,'FontWeight','bold');
%     ylabel('Counts','FontSize',16,'FontWeight','bold');
%     title(['Maximum hp position reached per trace, ' num2str(binsize) ' bp bins'], 'FontSize',16,'FontWeight','bold');
%     xlim([0 100]);
%    
%     h4=figure;
%     fig_handles = [fig_handles h4];
%     hold on;
%     for i=1:length(force)
%         binned_trace_max{i} = histc(mtracemax{i}, (0:binsize:100));
%         if isempty(binned_trace_max{i})
%             binned_trace_max{i} = zeros(1,size(binned_trace_max{i}, 1));
%         end
%         
%         legendstr{i} = [num2str((force(i)-forcebinw/2), '% .1f') ' to ' num2str(force(i)+forcebinw/2, '% .1f') 'pN ('  num2str(length(mtracemax{i}))  ' molecules)'];
%         plot((binsize/2:binsize:100+binsize/2), binned_trace_max{i}, 'Color', N(i, :));
%         
%     end
%     
%     xlabel('Max hairpin position reached (bp)','FontSize',16,'FontWeight','bold');
%     ylabel('Counts','FontSize',16,'FontWeight','bold');
%     title(['Maximum hp position reached per trace per force, ' num2str(binsize) ' bp bins'], 'FontSize',16,'FontWeight','bold');
%     legend(legendstr);
%     xlim([0 100]);
    
end

