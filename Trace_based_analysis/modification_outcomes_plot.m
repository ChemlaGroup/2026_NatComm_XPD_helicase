% From summary in bursts_classified.odt:


%% Block as separate category
% [HPconstruct,  block, advance+backslide, advance+rezip,
% advance+backslide/rezip combo, dissociation, clear]

burst_class = [9      4     9    5     1    2    35;...
               24     2     2     4    0    1    22;...
               15     5     6     8    1    2    30 ;...
               16     3    29     3    4    20     1 ;...
               14     7    15    20    3    6    36; ...
               17     1     3     2    1    4    15;...
               ];

           


%%
figure;
t = tiledlayout(2,3,'TileSpacing','compact');

% Create pie charts
ax1 = nexttile;
pie(ax1,burst_class(1, 3:end))
title(HPnames{9})

ax2 = nexttile;
pie(ax2,burst_class(2, 3:end))
title(HPnames{24})

ax6 = nexttile;
 pie(ax6,burst_class(3, 3:end))
 title(HPnames{15})

ax3 = nexttile;
pie(ax3,burst_class(4, 3:end))
title(HPnames{16})

ax4 = nexttile;
pie(ax4,burst_class(5, 3:end))
title(HPnames{14})


ax5 = nexttile;
pie(ax5,burst_class(6, 3:end))
title(HPnames{17})




% Create legend
%labels = {'block', 'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'clear'}
labels = {'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'clear'};
lgd = legend(labels);
lgd.Layout.Tile = 'south';

%%
figure; hold on;
bar(burst_class(:, 2:end), 'stacked');
labels = {'block', 'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'clear'};
lgd = legend(labels);

%% Bar plot, stacked, normalized

burst_class_norm = burst_class(:, 2:end)./sum(burst_class(:, 2:end), 2);

figure; hold on;
bar(burst_class_norm, 'stacked');
labels = {'block', 'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'clear'};
lgd = legend(labels);
ylabel('Fraction of bursts')
title('Fractions of bursts exhibiting different behaviours in encountering modification')

%%
burst_class_norm = burst_class(:, 3:end)./sum(burst_class(:, 3:end), 2);

figure; hold on;
bar(burst_class_norm, 'stacked');
labels = {'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'clear'};
lgd = legend(labels);
ylabel('Fraction of bursts')
title('Fractions of bursts exhibiting different behaviours in passing modification')


%% 'Block' bursts reclassified into other categories

% [HPconstruct,  backslide, rezip, backslide/rezip combo, dissociation, clear]
burst_class = [...
                9   9   8   1   3   35;...
                24  4   4   0   1   22;...
                16  32  3   4   20  1;...
                14  19  22  4   6   36;...
                17  3   3   1   4   15;...
                ];
            
            
%% Bar plot, stacked, normalized

burst_class_norm = burst_class(:, 2:end)./sum(burst_class(:, 2:end), 2);

figure; hold on;
barh([5 4 3 2 1],burst_class_norm, 'stacked'); % reversed y indeces to get top-down ordering
labels = {'backslide', 'rezip', 'backslide/rezip combo', 'dissociation', 'traversal'};
lgd = legend(labels);
title('Fractions of bursts exhibiting different behaviours in encountering modification')
xlabel('Fraction of bursts')

ylim([0.4 5.6])
const_labels{1} = num2str(burst_class(5, 1));
const_labels{2} = num2str(burst_class(4, 1));
const_labels{3} = num2str(burst_class(3, 1));
const_labels{4} = num2str(burst_class(2, 1));
const_labels{5} = num2str(burst_class(1, 1));

set(gca, 'YTickLabels', const_labels);

