% 200507, AVT -- major revision 250507

% Standalone function for calculating survival fraction, by burst
% -- Currently only implemented for unwinding.
% -- Can be used to calculate N-correction.

% 200521 Modified to calculate maximum of each burst; not rely on proc
% field, as it is not present in much of the data.

%221031: modified to display at step/2, in middle of bin.

%230206: modified to enable processivity cutoff, to effectively normalize survival probability by a minimum of interest. 
%241107: 
% -- modified to exlude bursts which begin above lower cutoff (e.g.
%after a reversal at higher hairpin position)
% -- implementing burst_list, subset of bursts to include (optional)

%250507: modifying to 
% a) correct binning error (shift in x)
% b) plot survival times directly, rather than histogramming.

%251205: modifying to round maximum and minimum to nearest basepair before
%thresholding

% Sample usage:
% [survival_frac, times, nbursts] = BurstSurvivalFrac(data{9}, 1, 30,[]);


function [survival_frac, drop_pos, nbursts]= BurstSurvivalFrac(data, makeplot, proc_cutoff, burst_list)

if nargin < 2
     makeplot = 0;
end

if nargin < 3
     proc_cutoff = 0;
end

if nargin<4
    burst_list = [];
end

max_proc = [];

[trace_num, burst_num]=BurstListConversion(data, burst_list);

%for i = 1:size(data,2)
    %for j=1:size(data{i}.burst2, 2)
 for i=trace_num
     
     if isempty(burst_num)
        Nburst = 1:size(data{i}.burst2, 2);
    else
        Nburst = burst_num{i}';
     end
    
     for j=Nburst
            %max_proc= [max_proc data{i}.proc];
            %entire trace
            %inds = find(data{i}.time>=data{i}.burst2(1, j) & data{i}.time<=data{i}.burst2(2, j));
            %unzipping only:
            inds = find(data{i}.time>=data{i}.burst2(1, j) & data{i}.time<=data{i}.burst2peak(j));
            burst_ext = data{i}.bp(inds);
            burst_max = max(burst_ext);
            burst_min = min(burst_ext);
            %burst_max = round(max(burst_ext));
            %burst_min = round(min(burst_ext));
            %if burst_min >proc_cutoff
            %    disp(['construct ' data{i}.construct  ', trace ' num2str(i) ', burst ' num2str(j) ', min ' num2str(burst_min)]);
            %end
            if round(burst_max) >= proc_cutoff && round(burst_min) <= proc_cutoff
                max_proc= [max_proc burst_max];
            end
        
    end
end

nbursts = length(max_proc);
drop_pos = [proc_cutoff sort(max_proc)];
survival_frac = (nbursts:-1:0)/nbursts;


if makeplot
   figure; 
   stairs(drop_pos, survival_frac);
   legendstr=[num2str(nbursts) ' bursts'];
   legend(legendstr);
   xlabel('Hairpin position');
   ylabel('Survival probability');
   if proc_cutoff
       title(['Survival probability for bursts reaching at least ' num2str(proc_cutoff) ' bp']);
   end
   
end


end