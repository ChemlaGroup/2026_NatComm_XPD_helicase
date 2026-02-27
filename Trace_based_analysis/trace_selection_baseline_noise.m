% 200122, AVT 

% Processing traces for selection by level of baseline noise.

% Locate the post-protein activity portion of the baseline in a trace, find
% the mean and standard deviation. Compare to the thresholds. If the trace does not pass the cutoff, try a
% constrained fit --  trace may still sufficiently well-behaved to be included
% based on constrained fit to histogram of extensions.

% Fit post-activity baseline to a Gaussian, within a constrained region,
% |extension| < 3bp by default.

%210603: the gaussian fit returns not std. dev, but sqrt(2)*std. dev. Not
% changing here for now, since my std. dev. cutoffs across the board would
% be affected and they are hard-wired into some of my codes. (The relative values are
% reasonable).

function [b3_mean, b3_stddev, b3fit_mean, b3fit_stddev] = trace_selection_baseline_noise(trace, baselinebinw, fit_constraint, b3_mean_threshold, b3_stddev_threshold)

if nargin<3
    fit_constraint =3;
    b3_mean_threshold = 2;
    b3_stddev_threshold = 1;
end

b3_ind = zeros(1, 2);

date = trace.date;
Fdiff = diff(trace.force);
Fchange_ind = find(abs(Fdiff(1:end-1)+Fdiff(2:end))>1.5)+1;

% Find interval corresponding to third baseline -- after XPD activity ends:
b3_ind(1) =find(abs(trace.time- trace.burst2(2, end)) < 0.5/trace.traprate);
% Find first timepoint with a large change in F after activity ceases,
% cut off an extra one to exlude the drop/rise.
temp_ind = Fchange_ind(Fchange_ind> b3_ind(1));

% If force not dropped down at end of trace, take the last point.
if isempty(temp_ind)
    temp_ind = length(trace.time);
end
b3_ind(2) = temp_ind(1)-1;
b3_ext = trace.bp(b3_ind(1):b3_ind(2));
b3_stddev = std(b3_ext);
b3_mean = mean(b3_ext);

% Possible that trace is still sufficiently well-behaved to be
% included based on constrained fit to histogram of extensions.
% Do fitting only if tether does not pass first cutoff.

if (abs(b3_mean) > b3_mean_threshold || b3_stddev > b3_stddev_threshold)
    
    baselinebins = floor(min(b3_ext)):baselinebinw:ceil(max(b3_ext));
    b3hist = histc(b3_ext, baselinebins);
    b3_fitpoints = b3hist(abs(baselinebins)<fit_constraint);
    baselinebins_fit = baselinebins(abs(baselinebins)<fit_constraint);
    
    try
        fit3=fit(baselinebins_fit'+baselinebinw/2, b3_fitpoints', 'gauss1');
        b3fit_stddev = fit3.c1;
        b3fit_mean = fit3.b1;
        %disp(['Fit: ' num2str(b3fit_mean), ', ' num2str(b3fit_stddev), ' ; element ' num2str(i)])
    catch
        disp(['Fitting |extension| < ' num2str(fit_constraint) ' failed for ' date '_' num2str(trace.refnums(4),'%03d') '. Attempting to fit full range.'])
        b3_fitpoints = b3hist;
        baselinebins_fit=baselinebins;
        try
            fit3=fit(baselinebins_fit'+baselinebinw/2, b3_fitpoints', 'gauss1');
            b3fit_stddev = fit3.c1;
            b3fit_mean = fit3.b1;
            %disp(['Unconstrained fit: ' num2str(b3fit_mean), ', ' num2str(b3fit_stddev), ' ; element ' num2str(i)])
        catch
            disp('Full range fit failed. Skipping fitting.')
            b3fit_stddev = [];
            b3fit_mean = [];
        end
        
    end
    
    % If first cutoff was passed, fit parameters are returned empty:
else
    b3fit_mean = [];
    b3fit_stddev = [];
end

end