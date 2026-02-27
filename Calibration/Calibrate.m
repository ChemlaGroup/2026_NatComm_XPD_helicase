% This is the shared driver code for doing calibrations using either the
% old trap or the fleezers. DO NOT MODIFY THIS CODE!
%
% This function returns a structure variable 'calData' which contains:
%     - bead: bead diameters for bead A and B
%     - alpha: the V to nm conversion factor for each bead in X and Y
%     - kappa: trap stiffness for each bead in X and Y
%     - offset: QPD offset for each bead in X and in Y
% The input variables are:
%     - Date: the data at which the data was acquired, e.g. '171217'
%     - rootfile: file identifier, e.g. '171217_002.dat'
%     - calarray: the parameters to get the calibration file, e.g. calparams = [900    880     125000  20000    8000    1200     100];
%     - showplots: an option to plot the calibration file, 0 - do not plot, 1 - plot
%     - writecal: an option to save the calibration file, 0 - no not save, 1 - save

% 181126: modified by AT to include Savefigdir as optional argument.
% 181126: for titles, replacing 'round' with num2str with precision
% specified -- my matlab version does not allow rounding to a precision.
% 221026: fixed bug I'd introduced at some point which separated generating
% plots and saving; if tried saving gcf without generating plot, ended up
% with something unrelated saved.
% 230201: undo writing cal, fit structures to base workspace if not on trap
% computer!
% 240208: adding TitleFontSizeMultiplier to figure titles, to prevent them
% overlapping (12 pt font size gets replaced by system default when figure
% reopened)

function [calData] = Calibrate(Date,rootfile,calarray,showplots,writecal, Savefigdir)
if nargin == 0
    Date = '240915'; % this may be changed
    datafnumber = 3; % this may be changed
    rootfile = [Date '_' num2str(datafnumber,'%03d') '.dat'];
end

global datadirectories          % parent directory for all data files
global trap_computer            % 1 if using instrument computer, 0 if not
global calparams
global fbslash

if nargin < 3
    prompts = {'Bead A diameter (nm)', 'Bead B diameter (nm)', 'f_{samp} (Hz)',...
        'f_{XY High} (Hz)', 'f_{Sum High} (Hz)', 'f_{Low} (Hz)', 'Averaging window'};
    defaults = {num2str(calparams(1)), num2str(calparams(2)), num2str(calparams(3)),...
        num2str(calparams(4)), num2str(calparams(5)), num2str(calparams(6)),...
        num2str(calparams(7))};
    calarray = getnumbers('Enter parameters:', prompts, defaults);
end
if nargin < 4
    showplots = 1;
end
if nargin < 5
    writecal = 1;
end

%AT, 181126
if nargin < 6
   Savefigdir =  [datadirectories Date fbslash]; % i.e. same as startpath
end

% initialize cal file name and data path
startpath = [datadirectories Date fbslash];
calfilename = ['cal' rootfile(1:(end-4))];
fitfilename = ['fit' rootfile(1:(end-4))];

% read the raw data
data = ReadJointFile(startpath,rootfile);

% Change array of input parameters to structure variable (better for saving
% user data later)
calpar.beadA = calarray(1);
calpar.beadB = calarray(2);
calpar.f_samp = calarray(3);
calpar.f_xyhi = calarray(4);
calpar.f_sumhi = calarray(5);
calpar.f_low = calarray(6);
calpar.avwin = calarray(7);

%which signal sets to analyze
calsets(1) = 1;%Trap A
calsets(2) = 1;%Trap B
calsets(3) = 0;%Detection laser

param = zeros(7,1);
param(7) = calpar.avwin; % ave windows
param(5) = calpar.f_sumhi; % f sum high

if calsets(1) == 1 %get cal for trap 1 (A)
    param(1) = calpar.beadA;
    param(4) = calpar.f_xyhi;% fxy high
    param(6) = calpar.f_low;% f low
    fitmodel = 'AliasedHydro';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.A_X;
    data_sub.Y = data.A_Y;
    data_sub.Sum = data.A_Sum;

    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    calData.beadA = param(1);
    calData.alphaAX = cal.alphaAX;
    calData.alphaAY = cal.alphaAY;
    calData.kappaAX = cal.kappaAX;
    calData.kappaAY = cal.kappaAY;
    calData.AXoffset = cal.AXoffset;
    calData.AYoffset = cal.AYoffset;

    allfit.fA = fit.f;
    allfit.AXSpecRaw = fit.AXSpecRaw;
    allfit.AYSpecRaw = fit.AYSpecRaw;
    allfit.predictedAX = fit.predictedAX;
    allfit.predictedAY = fit.predictedAY;
end

if calsets(2) == 1 %get cal for trap 2 (B)
    param(1) = calpar.beadB;
    param(4) = calpar.f_xyhi;% fxy high
    param(6) = calpar.f_low;% f low
    fitmodel = 'AliasedHydro';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.B_X;
    data_sub.Y = data.B_Y;
    data_sub.Sum = data.B_Sum;

    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    % All "cal" fields are named as if they were bead A from the
    % "CalibrateExact_Stripped" function, but they need to be renamed here,
    % since this is actually the calibration for bead B.
    calData.beadB = param(1);
    calData.alphaBX = cal.alphaAX;
    calData.alphaBY = cal.alphaAY;
    calData.kappaBX = cal.kappaAX;
    calData.kappaBY = cal.kappaAY;
    calData.BXoffset = cal.AXoffset;
    calData.BYoffset = cal.AYoffset;

    allfit.fB = fit.f;
    allfit.BXSpecRaw = fit.AXSpecRaw;
    allfit.BYSpecRaw = fit.AYSpecRaw;
    allfit.predictedBX = fit.predictedAX;
    allfit.predictedBY = fit.predictedAY;
end

if calsets(3) == 1 %get cal for detection laser on trap
    param(1) = 1000;%Bead B
    param(4) = 15e3;%fxy high
    param(6) = 200;%f low
    fitmodel = 'AliasedFiltered';
    data_sub.sampperiod = data.sampperiod;
    data_sub.X = data.C_X;
    data_sub.Y = data.C_Y;
    data_sub.Sum = data.C_Y_Sum;
    
    [cal,fit] = CalibrateExact_Stripped(param,fitmodel,startpath,rootfile,data_sub);
    
    calData.beadC = param(1);
    calData.alphaCX = cal.alphaAX;
    calData.alphaCY = cal.alphaAY;
    calData.kappaCX = cal.kappaAX;
    calData.kappaCY = cal.kappaAY;
    calData.CXoffset = cal.AXoffset;
    calData.CYoffset = cal.AYoffset;
    
    allfit.fC = fit.f;
    allfit.CXSpecRaw = fit.AXSpecRaw;
    allfit.CYSpecRaw = fit.AYSpecRaw;
    allfit.predictedCX = fit.predictedAX;
    allfit.predictedCY = fit.predictedAY;
    
end

if trap_computer
    assignin('base',calfilename,calData);
    assignin('base',fitfilename,allfit);
end

if (showplots || writecal) % 221026: if saving cal file, need to generate plots regardless of whether to leave them up or not.
    
    columns = sum(calsets);
    
    if ~trap_computer
        windowstyle = get(gcf, 'WindowStyle'); %AT, 220217; to avoid throwing warning when windowstyle is docked.
        if ~strcmp(windowstyle, 'docked')
            if sum(calsets) == 3
                figure('Position',[222          89        1626         871])
            else
                figure('Position',[25          49        900         600])
            end
        else
            figure;
        end
    else
        figure('Position',[80  100   1100   800])
    end
    
    
    set(gcf,'Name',['Calibration for ' rootfile(1:end-4)]);
    prev = 0;
    if calsets(1) == 1
        subplot(2,columns,1)
        loglog(allfit.fA, allfit.AXSpecRaw); hold on;
        loglog(allfit.fA, allfit.predictedAX, 'k');
        %title(['AX: \kappa = ' num2str(round(calData.kappaAX,3)) ', \alpha = ' ...
        %num2str(round(calData.alphaAX,2)) ', Offset = ' num2str(round(calData.AXoffset,2))]);
        title(['AX: \kappa = ' num2str(calData.kappaAX,3) ', \alpha = ' ...
            num2str(calData.alphaAX,4) ', Offset = ' num2str(calData.AXoffset,'%.3f')]);
        set(gca, 'TitleFontSizeMultiplier', 0.75);
        %, 'FontSize', 12); %, 'FontWeight', 'normal');
        % ylim([1e-9 5e-8])
        ylabel('Power (V^2 s)')
        subplot(2,columns,1+columns)
        loglog(allfit.fA, allfit.AYSpecRaw); hold on;
        loglog(allfit.fA, allfit.predictedAY, 'k');
        %title(['AY: \kappa = ' num2str(round(calData.kappaAY,2)) ', \alpha = ' ...
        %   num2str(round(calData.alphaAY,2)) ', Offset = ' num2str(round(calData.AYoffset,2))]);
        title(['AY: \kappa = ' num2str(calData.kappaAY,3) ', \alpha = ' ...
            num2str(calData.alphaAY,4) ', Offset = ' num2str(calData.AYoffset,'%.3f')]);
         set(gca, 'TitleFontSizeMultiplier', 0.75);
        %, 'FontSize', 12); %, 'FontWeight', 'normal');
        % ylim([1e-9 5e-8])
        prev = 1;
        xlabel('Frequency (Hz)')
        ylabel('Power (V^2 s)')
    end
    
    if calsets(2) == 1
        subplot(2,columns,prev+1)
        loglog(allfit.fB, allfit.BXSpecRaw); hold on;
        loglog(allfit.fB, allfit.predictedBX, 'k');
        %title(['BX: \kappa = ' num2str(round(calData.kappaBX,2)) ', \alpha = ' ...
        %   num2str(round(calData.alphaBX,2)) ', Offset = ' num2str(round(calData.BXoffset,2))]);
        title(['BX: \kappa = ' num2str(calData.kappaBX,3) ', \alpha = ' ...
            num2str(calData.alphaBX,4) ', Offset = ' num2str(calData.BXoffset,'%.3f')]);
         set(gca, 'TitleFontSizeMultiplier', 0.75);
        %, 'FontSize', 12); %, 'FontWeight', 'normal');
        % ylim([1e-9 5e-8])
        subplot(2,columns,prev+1+columns)
        loglog(allfit.fB, allfit.BYSpecRaw); hold on;
        loglog(allfit.fB, allfit.predictedBY, 'k');
        % ylim([1e-9 5e-8])
        %title(['BY: \kappa = ' num2str(round(calData.kappaBY,2)) ', \alpha = ' ...
        %   num2str(round(calData.alphaBY,2)) ', Offset = ' num2str(round(calData.BYoffset,2))]);
        title(['BY: \kappa = ' num2str(calData.kappaBY,3) ', \alpha = ' ...
            num2str(calData.alphaBY,4) ', Offset = ' num2str(calData.BYoffset,'%.3f')]);
         set(gca, 'TitleFontSizeMultiplier', 0.75);
        %, 'FontSize', 12); %, 'FontWeight', 'normal');
        prev = prev + 1;
        xlabel('Frequency (Hz)')
    end
    
    if calsets(3) == 1
        subplot(2,columns,prev+1)
        loglog(allfit.fC, allfit.CXSpecRaw); hold on;
        loglog(allfit.fC, allfit.predictedCX, 'k');
        %title(['CX: \kappa = ' num2str(round(calData.kappaCX,2)) ', \alpha = ' ...
        %num2str(round(calData.alphaCX,2)) ', Offset = ' num2str(round(calData.CXoffset,2))]);
        title(['CX: \kappa = ' num2str(calData.kappaCX,3) ', \alpha = ' ...
            num2str(calData.alphaCX,4) ', Offset = ' num2str(calData.CXoffset,'%.3f')]);
         set(gca, 'TitleFontSizeMultiplier', 0.75);
        %, 'FontSize', 12); %, 'FontWeight', 'normal');
        subplot(2,columns,prev+1+columns)
        loglog(allfit.fC, allfit.CYSpecRaw); hold on;
        loglog(allfit.fC, allfit.predictedCY, 'k');
        title(['CY: \kappa = ' num2str(calData.kappaCY,3) '. \alpha = ' ...
            num2str(calData.alphaCY,4) ', Offset = ' num2str(calData.CYoffset,'%.3f')]);
         set(gca, 'TitleFontSizeMultiplier', 0.75);
         %, 'FontSize', 12);%, 'FontWeight', 'normal');
    end
    
    %     ud.param = calpar;
    %     ud.fitted = calData;
    %     set(gcf, 'UserData', ud) % Save input parameters and fitted values into fig file
 
    if writecal
         saveas(gcf,[Savefigdir rootfile(1:(end-4)) '_cal.fig']);
    end
    
    if ~showplots
        close(gcf);
    end
    
end    
    
    if trap_computer && writecal % write calibration data for labview to grab later
    %if writecal
        %if showplots
        %saveas(gcf,[startpath rootfile(1:(end-4)) '_cal.fig']);
        saveas(gcf,[Savefigdir rootfile(1:(end-4)) '_cal.fig']);
        %end
        
        if data.oldtrap
            calout = [calData.alphaAX calData.alphaBX calData.kappaAX calData.kappaBX calData.AXoffset calData.BXoffset calData.alphaAY calData.alphaBY calData.kappaAY calData.kappaBY calData.AYoffset calData.BYoffset];
        else
            calout = [calData.alphaAX calData.alphaBX calData.kappaAX calData.kappaBX calData.AXoffset calData.BXoffset];
        end
        fid2 = fopen([startpath 'current.cal'],'w','ieee-be');
        fwrite(fid2,calout,'float64');
        fclose(fid2);
    end

end
