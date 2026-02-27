% 250718
% Takes in a cell array containing lists of start and end burst positions, and generates a survival curve.
% Useful for rezipping where bust classification and identification of bona fide rezipping needs to be done manually.

% burst_list{c}{i} has the form:
% date trace_number burst_number start_position end_position

% BurstSurvivalFracRezipping_plotDatasets(rz_survival_list_3set, [19 22], 55);

function BurstSurvivalFracRezipping_plotDatasets(burst_list, constructs, proc_cutoff)

global HPcolours
global HPnames

if nargin < 2
    error('Please select constructs.')
end
    
if nargin < 3
     proc_cutoff = 55;
end

figure; hold on;

for c=constructs
    [survival_frac{c}, drop_pos{c}, nbursts{c}]= BurstSurvivalFracRezipping(burst_list{c}, proc_cutoff, 0);
    stairs(drop_pos{c}, (survival_frac{c}*nbursts{c}), 'Color', HPcolours{c}, 'LineWidth', 2, ...
     'DisplayName', [HPnames{c} ', ' num2str(nbursts{c}) ' bursts' ]);
end

 xlabel('Hairpin position');
   %ylabel('Survival probability');
   ylabel('Surviving bursts');
   if proc_cutoff
       %title(['Survival probability for bursts rezipping from ' num2str(proc_cutoff) ' bp']);
       title(['Survival for bursts rezipping from ' num2str(proc_cutoff) ' bp']);
   end
   xlim([0 proc_cutoff]);
   set ( gca, 'xdir', 'reverse' )
   %ylim([0 1]);
   
end


%[survival_frac, drop_pos]= BurstSurvivalFracRezipping(hp19rz, 55,1);
function [survival_frac, drop_pos, nbursts]= BurstSurvivalFracRezipping(burst_list, proc_cutoff, makeplot)

if nargin < 2
    proc_cutoff = 55;
end

if nargin < 3
    makeplot = 0;
end


%[trace_num, burst_num]=BurstListConversion(data, burst_list(:,1:3));
end_proc = [];

for i= 1:size(burst_list,1)
    if burst_list(i,4) >= proc_cutoff && burst_list(i,5) <= proc_cutoff
        end_proc = [end_proc burst_list(i,5)];
        
    end
end

nbursts = length(end_proc);

drop_pos = [proc_cutoff sort(end_proc, 'descend')];
survival_frac = (nbursts:-1:0)/nbursts;

if makeplot
    figure;
    stairs(drop_pos, survival_frac);
    legendstr=[num2str(nbursts) ' bursts'];
    legend(legendstr);
    xlabel('Hairpin position');
    ylabel('Survival probability');
    if proc_cutoff
        title(['Survival probability for bursts rezipping from ' num2str(proc_cutoff) ' bp']);
    end
    xlim([20 proc_cutoff]);
    set ( gca, 'xdir', 'reverse' )
    ylim([0 1]);
end



end