% AVT
% processing constant force files, analogous to my CF_analysis_raw code,
% but using Monika's functions.

% First-pass analysis; contains some ad-hoc criteria to estimate baseline,
% etc. In particular, extension calculations for whole trace are based on
% force within a window at the start: will not be accurate if force is raised afterwards.
% Note also that since the estimated offset is DNA extension at
% 'measurement' force, extension at low force will be very inaccurate with
% offset subtraction.

%220722: modified a==b to abs(a-b)<eps*10 structure for direction
%determination. Came across a file where old structure failed.

%220617: modified lower ad-hoc force cutoff to 7 pN. My low force in some
%files had gone over 5 pN and this was throwing off the calculations.
% Added makecalplot parameter. If plotting force-extension curves
% automatically, can avoid plotting same calibration again.

%230725: modifying to take forces in both directions into account.
% Especially important given diagonal pulling, from misalignment of
% steerable mirror (accidental tilt).

%230918: if first baseline absent, subtract median extension over ~right force range as baseline estimate.

%[DNALengthbp, FM] = CF_plot(Date, calfilenumber, offsetfilenumber, CFfilenumber, 0, 0)
% [~, ~]= CF_plot('250107', 36, 37, 49, 1, 1);

function  [DNALengthbp, FM] = CF_plot(Date, calfilenumber, offsetfilenumber, CFfilenumber, makeplot, makecalplot)

global datadirectories
global fbslash
global calparams
global construct
global analysis

cd(analysis);

Savefigdir = [datadirectories Date '_2Doffset' fbslash];
if ~exist(Savefigdir, 'dir')
    mkdir(Savefigdir);    
end

% for ease of sorting later
if ~exist([Savefigdir 'Not_dipped' fbslash], 'dir')
    mkdir([Savefigdir 'Unusable_or_noactivity_tethers' fbslash]);
    mkdir([Savefigdir 'Not_dipped' fbslash]);
end

% if nargin <5
%     usetrapF = 'AB'; % by default, use forces from both traps
% end

 if nargin <5
     makeplot = 1;
 end
 
 if nargin <6
 makecalplot = 1; % 1 - plot, 0 - do not plot
 end
 
plotFluo=0;
 pullingCurve=0;
 
%construct = 'Hairpin-LambdaRH-dT-10-86bp';


param.Date = Date;
param.calfilenumber = calfilenumber;
param.CFnumber = CFfilenumber;

param.offset_file = 1; % 1 - use offset file, 0 - subtract only QPD offset
param.file = [param.Date fbslash '_' num2str(param.CFnumber,'%03d')];

startpath = [datadirectories param.Date fbslash];

% calibrate
param.plotCal = makecalplot; % 1 - plot, 0 - do not plot
calRootFile = [param.Date '_' num2str(param.calfilenumber,'%03d') '.dat'];
calfilename = ['cal' calRootFile(1:(end-4))];
calData = Calibrate(param.Date,calRootFile,calparams,param.plotCal, 1, Savefigdir);

% get raw data
CFRootFile = [param.Date '_' num2str(param.CFnumber,'%03d') '.dat'];
CFfilename = ['data_' CFRootFile(1:(end-4))];
CFData = ReadJointFile(startpath,CFRootFile);
param.instrument = CFData.oldtrap ; % 1 - old trap, 0 - fleezers
%if CFData.scandir == 0
if abs(CFData.t1x-CFData.t2x) < eps*10   %CFData.t1x == CFData.t2x
    param.direction = 'Y'; % Y or X
    %elseif CFData.scandir == 1
elseif abs(CFData.t1y-CFData.t2y) < eps*10                                %CFData.t1y == CFData.t2y
    param.direction = 'X'; % Y or X
end

disp([num2str(CFfilenumber) ': ' param.direction]);

if param.instrument
    % get fixedxy results from bead grid analysis -- for old trap
    beta = getBetaFactors(startpath, param.Date); disp('Using beta from file.')
else
    beta = eye(2); % use dummy identity beta for fleezers
end

CFData = trap_data_nmConversion(CFData, param, beta); % convert from V or MHz to nm
% get offset data
% offsetfilenumber = param.calfilenumber + 1;
offsetRootFile = [param.Date '_' num2str(offsetfilenumber,'%03d') '.dat'];
offsetfilename = ['offset' offsetRootFile(1:(end-4))];
offsetData = ReadJointFile(startpath, offsetRootFile);
offsetData = trap_data_nmConversion(offsetData, param, beta); % convert from V or MHz to nm

% use 2D offset subtraction:
trapData = force_offset_subtraction(offsetData, CFData, calData, param); % using XY offset

% % Taking both x and y into account
FA = sqrt(trapData.force_AX.^2 + trapData.force_AY.^2);
FB = sqrt(trapData.force_BX.^2 + trapData.force_BY.^2);
FM = (FA + FB)./2;
DNAext=trapData.DNAext;
 
% FX = (trapData.force_AX+trapData.force_BX)./2;
% FY = (trapData.force_AY+trapData.force_BY)./2;
% FM = sqrt(FX.^2+FY.^2); %Ftot in my notes
% DNAext = sqrt(trapData.DNAext_X.^2 + trapData.DNAext_Y.^2);



% Find the mean force value from the meaningful part of the file -- exclude
% e.g. forces after tether rupture and before/after the force is turned up.
% A bit ad-hoc and somewhat arbitrary cutoffs here.
% Add time constraint, since I sometimes step force down later in the
% trace.
relevant_inds = find(FM > 7 & FM<15 &CFData.time<60);
mean_force = mean(FM(relevant_inds));
median_force = median(FM(relevant_inds));

[mean_force, median_force];

%% load construct parameters
ConstructConst = ConstructConstants(construct);

hp_params(1)= ConstructConst.dsDNA1; % 3.055e3; %length of ds handles
hp_params(2) =ConstructConst.ssDNA1 ; %length of ss loading site
hp_params(3) = 86; %89; %length of hairpin
hp_params(4) = 4; %length of T loop in hairpin

%instantaneous force:
DNALengthbp = bpunwound(FM,DNAext, mean_force,hp_params);
% average force:
%DNALengthbp = bpunwound(mean_force*ones(1, numel(DNAext)),DNAext, mean_force,hp_params);

%% Plots -- need to be converted to correct variable names

% Select color by force
forces = [6:1:14];

N = [0.8 1 0.4;...
    1 0.8 0;...
    1 0 0 ;...
    1 0.6 0.6;
    1 0 0.8;...
    0.6 0 0.6 ; ...
    0 0 1;...
    0 0.6 0.8; ...
    0 0.8 1];

[stgTime, stgPos] = StageMoveTime(Date, [param.Date '_' num2str(param.CFnumber,'%03d') '_stagetime.txt']);
force_colour_ind = find(forces == round(mean_force));
if isempty(force_colour_ind)
   force_colour = [1, 1, 1];
else
   force_colour= N(force_colour_ind, :);
end
% ad-hoc baseline correction:
%ext_baseline_estimate = median(DNALengthbp(relevant_inds));

% Baseline: median extension before moving out of lower channel
ext_baseline_estimate = median(DNALengthbp(CFData.time < stgTime(2)));

if isnan(ext_baseline_estimate)
    % estimate based on overall correct force range -- in case first
    % baseline is missing
    ext_baseline_estimate= median(DNALengthbp(FM > 7 & FM<15));
end

DNALengthbp = DNALengthbp-ext_baseline_estimate;

if makeplot
    % Plotting extension
    figure;
    set(gcf, 'FileName',[Savefigdir Date '_' num2str(param.CFnumber,'%03d')],...
        'Name', [Date '_' num2str(param.CFnumber, '%03d') '_forces_length']);
    
    length = subplot(4,1, [1 2]); hold on;
    %plot(CFData.time, DNALengthbp, 'Color', N(force_colour_ind, :));
    plot(CFData.time, DNALengthbp, 'Color', force_colour);
    legend([Date '_' num2str(param.CFnumber, '%03d') '; ' num2str(mean_force, '%0.1f') ' pN']);
    %annotation('textbox', [0.15, 0.83, 0.162, 0.085], 'String',[num2str(ext_baseline_estimate) 'bp baseline correction'],'FitBoxToText','on');
    legend('off')
    ylabel('Length (bp)');
    %xlabel('Time (s)'); %buried by lower plot
    
    
    %ylim([ext_baseline_estimate-10  ext_baseline_estimate+100]);
    
    %Default setting:
    ylim([-10 100]);
    
    %ylim([-10 30]);
    
    % Time in upper stream: between third and fourth entry
    % Moves start and end at 2nd and 5th
    
    [~, t2_ind] =min(abs(CFData.time - stgTime(2)));
    
    if size(stgTime, 1)>=4
        [~, t4_ind] =min(abs(CFData.time - stgTime(4)));
    else
        [~, t4_ind] = max(CFData.time);
    end
        
    if size(stgTime, 1)>=5
        [~, t5_ind] =min(abs(CFData.time - stgTime(5)));
    else
        [~, t5_ind] = max(CFData.time);
    end
    %[~, t3_ind] =min(abs(CFData.time - stgTime(3)));
    %[~, t4_ind] =min(abs(CFData.time - stgTime(4)));
    
    plot([CFData.time(t2_ind) CFData.time(t5_ind)],[DNALengthbp(t2_ind), DNALengthbp(t5_ind)], 'go', 'LineStyle', 'none');
    plot(CFData.time(t4_ind), DNALengthbp(t4_ind), 'ko', 'LineStyle', 'none');
    %plot(CFData.time(t5_ind),DNALengthbp(t5_ind), 'ro');
    
    
    forces = subplot(4,1,3, 'align');
    plot(CFData.time, FA, 'r', CFData.time, FB, 'g', CFData.time, FM, 'b');
    %plot(CFData.time, FA, 'r', CFData.time, FB, 'g', CFData.time, FM_AB, 'b'); hold on;
    %plot(CFData.time, FX, 'c', CFData.time, FY, 'm', CFData.time, FM, 'k');
    ylabel('Force (pN)');
    %xlabel('Time (s)');
    
    %legend('Trap A', 'Trap B', 'Average (A and B)', 'FX', 'FY', 'Total (X and Y)');
    
    legend('Trap A', 'Trap B', 'Average (A and B)');
    ylim([5 15]);
    
    plot3 = subplot(4,1,4);
    
    if plotFluo
        plot(CFData.apdtime, CFData.apd1, 'g'); set(g,'DisplayName','Cy3'); hold on;
        plot(downsample(CFData.apdtime,10), downsample(CFData.apd1,10), 'Color', [0,0.75,0]);
        hold on;
        plot(CFData.apdtime, data.apd2, 'r'); set(r,'DisplayName','Cy5'); hold on;
        plot(downsample(CFData.apdtime,10), downsample(CFData.apd2,10),'Color',[0.45,0,0]);
        ylabel('Photon Rate (kHz)')
        %axis tight
    else        %Plot bead positions in orthogonal direction -- to determine when move takes place
        
        if param.direction == 'X'
            plot(CFData.time, CFData.A_Y, 'r');
            hold on;
            plot(CFData.time, CFData.B_Y, 'b');
            legend('A_Y', 'B_Y');
            ylabel(['Bead Position in Y (Hz)'])
        elseif param.direction == 'Y'
            plot(CFData.time, CFData.A_X, 'r');
            hold on;
            plot(CFData.time, CFData.B_X, 'b');
            legend('A_X', 'B_X');
            ylabel(['Bead Position in X (Hz)'])
        end
        
           xlabel('Time (s)')
        
    end
    
    linkaxes([length, forces, plot3], 'x')
    
    %saveas(gcf, [startpath Date '_' num2str(param.CFnumber,'%03d') '_forces_A1.fig'])
    saveas(gcf, [Savefigdir Date '_' num2str(param.CFnumber,'%03d') '_forces_A1.fig'])
    
    
    if pullingCurve
        figure;
        set(gcf, 'Name', [Date '_' num2str(param.CFnumber, '%03d') '_force_extension']);
        plot(DNAext, FA, 'r', DNAext, FB, 'g', DNAext, FM, 'b')
        ylabel('Force (pN)');
        xlabel('Extension (nm)');
        legend('Trap A', 'Trap B', 'Average');
        hold on;
        axis([800, 1300, 0, 25])
        axis 'auto y'
        %saveas(gcf, [startpath Date '_' num2str(param.CFnumber,'%03d') '_force_extension_A1.fig'])
        saveas(gcf, [Savefigdir Date '_' num2str(param.CFnumber,'%03d') '_force_extension_A1.fig'])
    end

end

end