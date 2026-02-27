
% Front end to CF_plot

% To generate model first:
% [HPF_av,HPxtot,hpGtot] = ModelHairpinUnwind(9,'TTTT',10);hpmodel = [HPxtot; HPF_av];

%plot FECs together
%[fake_offsets_all, ~] = FEC_overlay_processing(flist,1, 'XPD buffer');xlim([850 1200]); ylim([0 20]);

for i=1:size(flist, 1)
   %Plot FECs individually
%    if exist('hpmodel', 'var')        
%        plot_force_extension(num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), [5 10], 'AB', hpmodel);
%    else
%         disp('No unfolding model pre-generated!')
%        plot_force_extension(num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), [5 10]);
% end
    if ~isempty(flist(i, 5)) && ~isnan(flist(i, 5))
        [~, ~]=CF_plot(num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 5),1, 1);
        %don't plot calibrations twice
    end
end

