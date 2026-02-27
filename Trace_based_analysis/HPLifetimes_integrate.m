% 230212, AVT

% Sum over the mean lifetimes calculated by HPlifetimes_errors over a particular
% range

% 250627: 

% Sample usage
% [m_hist, sem, position_hists, bins, avgforces]=HPlifetimes_errors(0, 1, 1, 25, 'Mean', 1,data);
% [lifetime_sums, cumsum_errors]= HPLifetimes_integrate(m_hist, sem, bin_centers, [18 37], [9 24 15 14 16 17]);

function [lifetime_sums, cumsum_errors] = HPLifetimes_integrate(lifetime_hist, errors, bin_centers, range, constructs)

global HPcolours
global HPnames

barplot=1;

bin_inds(1) = find(bin_centers == range(1));
bin_inds(2) = find(bin_centers == range(2));
labels = {};

nbins = bin_inds(2)-bin_inds(1)+1;
lifetime_sums = NaN(numel(constructs), nbins);
cumsum_errors = lifetime_sums;


for i=1:numel(constructs)
    lifetime_sums(i, :)=cumsum(lifetime_hist{constructs(i)}(bin_inds(1):(bin_inds(2))));
    % earlier version, when cumulative time over single range calculated:
    %lifetime_sums(i)=sum(lifetime_hist{constructs(i)}(bin_inds(1):bin_inds(2)));
    % adding in quadrature:
    cumsum_errors(i, :) = sqrt(cumsum(errors{constructs(i)}(bin_inds(1):(bin_inds(2))).^2));
    % upper limit of error -- adding errors linearly:
    %cumsum_errors(i, :) = cumsum(errors{constructs(i)}(bin_inds(1):(bin_inds(2))));
    labels = [labels HPnames{constructs(i)}];
end

figure; hold on;

%X=categorical(labels);
%X = reordercats(X, labels);
%bar(X, lifetime_sums);

xvals = (bin_centers(bin_inds(1)):bin_centers(bin_inds(2))); %+(bin_centers(2)-bin_centers(1))/2;

for i=1:numel(constructs)
    plot(xvals, lifetime_sums(i, :), 'Color', HPcolours{constructs(i)});
end

for i=1:numel(constructs)
    errorbar(xvals,lifetime_sums(i, :), cumsum_errors(i, :), 'LineStyle', 'none', 'Marker', 'o', 'Color', HPcolours{constructs(i)});
end

legend([labels labels], 'Location', 'northwest');
xlim(range);
 %title(['Time spent over positions ' num2str(range) ' bp']);

 title(['Cumulative time spent between positions ' num2str(range(1)) ' and ' num2str(range(2)) '  bp']);

 ylabel('Cumulative time (s)')
 xlabel('Hairpin position (bp)')
 
 if barplot
     for c=1:numel(constructs)
         col(c, :)=HPcolours{constructs(c)};
     end
     
     figure; hold on;
     X=categorical(labels);
     X = reordercats(X, labels);
     b=bar(X, lifetime_sums(:, end));
     b.FaceColor = 'flat';
     b.CData = col;
     errorbar(X, lifetime_sums(:, end),cumsum_errors(:, end), 'k', 'LineStyle', 'none', 'LineWidth', 3, 'CapSize', 15);
     ylabel('Cumulative time (s)');
     title(['Cumulative mean time spent between positions ' num2str(range(1)) ' and ' num2str(range(2)) '  bp']);
 end
 
 
end