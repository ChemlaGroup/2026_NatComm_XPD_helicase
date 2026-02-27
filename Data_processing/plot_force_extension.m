% function plot_force_extension(trapData, construct, param)

% AT, 181126: introduced Savefigdir variable; save figures in separate directory; pnc - 'processed by new code'
% Calibrations also saved in this directory.

%AT, 181206: usetrapF parameter; in case one trap has something wrong with
%it (e.g. bad calibration), can use force from the other trap only instead of average force. Parameter should
%be 'A' or 'B' or 'AB'. Changed to an input on 190313.
% AT, 181210: manfitrange is an optional parameter to manually set the fitting range
% for the current molecule; should be of the form [4 12]).

% AT, 190313: plot hairpin model, optionally (requires it to be in memory
% already); pass in hpmodel, as [HPxtot; HPF_av]. Modified legend workings
% to optionally have another entry.

% AT, 230801: subtract offsets and calculate forces using both dimensions, not just dominant one.
% AT, 230813: updating to include improvements of code on the Old Trap.

% Example usage:
% plot_force_extension('181129',152, 153, 162)

function plot_force_extension(Date,calfilenumber, offsetfilenumber, fextnumber, manfitrange, usetrapF, hpmodel)
global datadirectories
global fbslash
global calparams
global construct
global WLC_param

Savefigdir = [datadirectories Date '_2Doffset' fbslash];
if ~exist(Savefigdir, 'dir')
    mkdir(Savefigdir);
end

if nargin <5
    manfitrange = [];
    usetrapF = 'AB'; % by default, use forces from both traps
    hpmodel =0;
end

if nargin <6
    usetrapF = 'AB'; % by default, use forces from both traps
    hpmodel =0;
end

if nargin <7
    hpmodel =0;
end


%% initialize if no input
if nargin == 0
    Date = '180616';
    calfilenumber =  266;
    fextnumber =    268;
    manfitrange = [];
    
    prompts = {'Construct', 'dsDNA interphosphate distance (nm/bp)', 'ssDNA interphosphate distance (nm/bp))', 'dsDNA persistance length (nm)',...
        'ssDNA persistance length (nm)', 'dsDNA stretch modulus (pN)', 'ssDNA stretch modulus (pN)', 'kT'};
    defaults = {construct, num2str(WLC_param.hds), num2str(WLC_param.hss),...
        num2str(WLC_param.Pds), num2str(WLC_param.Pss), num2str(WLC_param.Sss),...
        num2str(WLC_param.Sds), num2str(WLC_param.kT)};
%     calarray = getnumbers('Enter parameters:', prompts, defaults);
    
    fext_array = getnumbers_fext('Enter parameters:', prompts, defaults);
    construct = cell2mat(fext_array(1));
    WLC_param.hds = str2num(cell2mat(fext_array(2))); % [nm/bp] interphosphate distance 
    WLC_param.hss = str2num(cell2mat(fext_array(3))); % nm/nt
    WLC_param.Pds = str2num(cell2mat(fext_array(4))); % nm
    WLC_param.Pss = str2num(cell2mat(fext_array(5))); % nm
    WLC_param.Sss = str2num(cell2mat(fext_array(6))); % pN
    WLC_param.Sds = str2num(cell2mat(fext_array(7))); % pN
    WLC_param.kT = str2num(cell2mat(fext_array(8)));% Zhi's values (dsDNA)
    
end
%%
param.Date = Date;
param.calfilenumber = calfilenumber;
param.fextnumber = fextnumber;

param.offset_file = 1; % 1 - use offset file, 0 - subtract only QPD offset
param.file = [param.Date fbslash '_' num2str(param.fextnumber,'%03d')];

startpath = [datadirectories param.Date fbslash];

% calibrate
param.plotCal = 1; % 1 - plot, 0 - do not plot
calRootFile = [param.Date '_' num2str(param.calfilenumber,'%03d') '.dat']; 
calfilename = ['cal' calRootFile(1:(end-4))];
calData = Calibrate(param.Date,calRootFile,calparams,param.plotCal, 1, Savefigdir);



% get raw data
fextRootFile = [param.Date '_' num2str(param.fextnumber,'%03d') '.dat']; 
fextfilename = ['data_' fextRootFile(1:(end-4))];
fextData = ReadJointFile(startpath,fextRootFile);
param.instrument = fextData.oldtrap ; % 1 - Old Trap, 0 - fleezers

if param.instrument
    % get fixedxy results from bead grid analysis -- for old trap
    beta = getBetaFactors(startpath, param.Date); disp('Using beta from file.')
else
    beta = eye(2); %disp('Using identity beta') % use dummy identity beta for fleezers
end

if fextData.scandir == 0
param.direction = 'Y'; % Y or X
elseif fextData.scandir == 1
param.direction = 'X'; % Y or X
end
fextData = trap_data_nmConversion(fextData, param, beta); % convert from V or MHz to nm
% get offset data
% offsetfilenumber = param.calfilenumber + 1;
offsetRootFile = [param.Date '_' num2str(offsetfilenumber,'%03d') '.dat']; 
offsetfilename = ['offset' offsetRootFile(1:(end-4))];
offsetData = ReadJointFile(startpath, offsetRootFile);  
offsetData = trap_data_nmConversion(offsetData, param, beta); % convert from V or MHz to nm
% subtract offset
%param.offset_file = 0; %testing only
trapData = force_offset_subtraction(offsetData, fextData, calData, param); % using XY offset subtraction code


%% load construct parameters

ConstructConst = ConstructConstants(construct);

if ~isempty(manfitrange)
    ConstructConst.fitrange = [manfitrange];
end



%% plot the WLC model
set(0,'defaultfigureposition',[100  200  1050  659]') % Put figures in middle of laptop screen.
max_force = max([trapData.force_AX, trapData.force_AY, trapData.force_BX, trapData.force_BY]);
F = 0:0.05:max_force;
% [WLC_param] = WLC_parameters;
bp = 1*WLC_param.hds*XWLCContour(F,WLC_param.Pds,WLC_param.Sds,WLC_param.kT); % nm/bp
nt = 1*WLC_param.hss*XWLCContour(F,WLC_param.Pss,WLC_param.Sss,WLC_param.kT); % nm/nt
                
if ConstructConst.hairpin == 1
    % closed model
    WLC_ext = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp;
    % open model
    WLC_ext_op = ConstructConst.ssDNA2.*nt + ConstructConst.dsDNA2.*bp;
    figure
    plot(WLC_ext, F, '--r', 'LineWidth', 2, 'DisplayName', 'WLC closed')
    hold on
    plot(WLC_ext_op, F, '--k', 'LineWidth', 2, 'DisplayName', 'WLC open')
else
    WLC_ext = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp;
    figure
    plot(WLC_ext, F, '--r', 'LineWidth', 2, 'DisplayName', 'WLC model')
end

set(gcf, 'Name', ['FEC for' fextRootFile]);

%% plot the data
% Filter TrapSep
TrapSepFilt  = movmean(trapData.TrapSep,30);
%TrapSepFilt  = movingmean(trapData.TrapSep, 31, 2, 1);   % AT: input of this function for window size should be odd, so set 31
% remove jumps in trap position
TrapSepFilt_diff = diff(TrapSepFilt);
normalized_diff = diff(TrapSepFilt)/mean(TrapSepFilt_diff(TrapSepFilt_diff>0));
TrapSepFiltClean = TrapSepFilt(normalized_diff > 0.8 | normalized_diff < -0.8);
trapData.force_AX = trapData.force_AX(normalized_diff > 0.8 | normalized_diff < -0.8);
trapData.force_BX = trapData.force_BX(normalized_diff > 0.8 | normalized_diff < -0.8);
trapData.force_AY = trapData.force_AY(normalized_diff > 0.8 | normalized_diff < -0.8);
trapData.force_BY = trapData.force_BY(normalized_diff > 0.8 | normalized_diff < -0.8);
trapData.DNAext = trapData.DNAext(normalized_diff > 0.8 | normalized_diff < -0.8);
% find the top peaks
[numMAX, indMAX] = findpeaks(TrapSepFiltClean);
% find the bottom peaks
TrapSepInverse = [-max(numMAX) -TrapSepFiltClean -max(numMAX)];
[numMIN, indMIN] = findpeaks(TrapSepInverse);
indMIN = indMIN-1;
% count number of raster scans
%rasterNUM = (max([length(indMIN) length(indMAX)]))/2;
rasterNUM = length(indMAX);

FA = sqrt(trapData.force_AX.^2 + trapData.force_AY.^2);
FB = sqrt(trapData.force_BX.^2 + trapData.force_BY.^2);
FM = (FA + FB)./2;
DNAext=trapData.DNAext;

%Default: average force, usetrapF = 'AB'. Otherwise:
if usetrapF == 'A'
    FM = FA; disp('Using trap A force.');
elseif usetrapF == 'B'
    FM = FB; disp('Using trap B force.');
end

DNAext_offset = [];
% get fake offset
DNAext_offset = mean(DNAext(FM > ConstructConst.fitrange(1) & FM < ConstructConst.fitrange(2)));
WLC_offset = mean(WLC_ext(F > ConstructConst.fitrange(1) & F < ConstructConst.fitrange(2)));
fake_offset_0 = DNAext_offset - WLC_offset;

% SY modified code. Will comment above until SM meeting 012119 
fake_offset = getFakeoffset(DNAext, FM, ConstructConst.fitrange, fake_offset_0, construct) ;

% substract the fake offset
DNAext = DNAext-fake_offset;

for ii = 1:1:rasterNUM
    indxS = indMIN(ii);
    indxE = indMAX(ii);
    indxEE = indMIN(ii+1);
hold on
plot(DNAext(indxS:indxE),FM(indxS:indxE),'Color', [0.9100    0.4100    0.1700], 'LineWidth', 1.5, 'DisplayName', 'mean F pull')
hold on
plot(DNAext(indxE:indxEE),FM(indxE:indxEE), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5,'DisplayName', 'mean F relax')
hold on
plot(DNAext(indxS:indxE),FA(indxS:indxE),'-.','Color', [0.9100    0.4100    0.1700], 'LineWidth', 0.5,'DisplayName', 'FA pull')
hold on
plot(DNAext(indxE:indxEE),FA(indxE:indxEE),'-.', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5,'DisplayName', 'FA relax')
hold on
plot(DNAext(indxS:indxE),FB(indxS:indxE),'-.','Color', [0.9100    0.4100    0.1700], 'LineWidth', 0.5,'DisplayName', 'FB pull')
hold on
plot(DNAext(indxE:indxEE),FB(indxE:indxEE),'-.', 'Color', [0 0.4470 0.7410], 'LineWidth',0.5,'DisplayName', 'FB relax')
hold on
end

if hpmodel
   plot(hpmodel(1, :), hpmodel(2, :), 'b', 'LineWidth', 1.5,'DisplayName', 'HP unfolding model') 
end

xlim([800 1200])
ylim([0 20])
%xlim([min(DNAext)-100 max(DNAext)+100])
%ylim([min([0 min(FM)-5]) max(FM)+5])
%ylim([0 max(FM)+5])
xlabel('DNA extension (nm)')
ylabel('Force (pN)')
if param.instrument == 1 
    instrument = 'Old Trap';
elseif param.instrument == 0
    instrument = 'Fleezers';
end
if rasterNUM == 1
  legend('location','northwest')
end
legend off
title({[ 'Force vs. Extension, ' num2str(fake_offset) ' nm offset ' '(' param.Date '\_' num2str(param.fextnumber,'%03d') '), bead diameters (' num2str(calparams(1)) ', ' num2str(calparams(2)) ') nm']...
    [construct ' construct, ' instrument ', force from trap(s) ' usetrapF ' for fits in range (' num2str(ConstructConst.fitrange(1)) ', ' num2str(ConstructConst.fitrange(2)) ') pN']}, 'FontSize', 12)
% box on

% saveas(gcf, [startpath param.Date '_' num2str(param.fextnumber,'%03d') '_fext2.fig'])
saveas(gcf, [Savefigdir param.Date '_' num2str(param.fextnumber,'%03d') '_fext2.fig'])

end
