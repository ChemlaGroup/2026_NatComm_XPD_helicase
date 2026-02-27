
%231113, AVT
% Based on earlier code snippet, datasets_survival_probability_plots.

function BurstSurvivalFrac_plotDatasets(data, constructs, burst_list)

global HPcolours
global HPnames

%N = distinguishable_colors(6);

figure; hold on;

if nargin <3 || isempty(burst_list)
   burst_list =cell(1, max(constructs)); % empty placeholder --> all bursts get included
   %burst_list = []; 
end

for c=constructs
 [survival_frac{c}, drop_pos{c}, nbursts{c}]= BurstSurvivalFrac(data{c}, 0, 30, burst_list{c});
 stairs(drop_pos{c}, (survival_frac{c}*nbursts{c}), 'Color', HPcolours{c}, 'LineWidth', 2, ...
     'DisplayName', [HPnames{c} ', ' num2str(nbursts{c}) ' bursts' ]);
end

legend;
xlim([30 75]);
%set(gca, 'XGrid', 'on')
%set(gca, 'XMinorGrid', 'on')
%set(gca, 'YGrid', 'on')
%set(gca, 'YMinorGrid', 'on')

xlabel('Hairpin position'); ylabel('Surviving number of bursts')

% With processivity cutoff (to get normalization by position)

figure; hold on;

for c = constructs
    stairs(drop_pos{c}, survival_frac{c}, 'Color', HPcolours{c}, 'LineWidth', 2, ...
        'DisplayName', [HPnames{c} ', ' num2str(nbursts{c}) ' bursts' ]);
end

plot([35 35], [0 1], '--k','DisplayName', 'Modification site');
legend;
xlim([30 75]);
ylim([0 1]);
% set(gca, 'XGrid', 'on')
% set(gca, 'XMinorGrid', 'on')
% set(gca, 'YGrid', 'on')
% set(gca, 'YMinorGrid', 'on')
%f=get(gca, 'Children');
%legend(f(end:2), legendstr); %doesn't work with DisplayName setup that I
%can see for now
%xlabel('Hairpin position (bp)', 'FontWeight', 'bold'); ylabel('Surviving fraction', 'FontWeight', 'bold')
xlabel('Hairpin position (bp)'); ylabel('Surviving fraction')

end