% AVT

% This code version takes into account X and Y offsets.
% Modifying plot_force_extension code to return DNA extension and force for
% further analysis.

% 230917: cleaning up FA and FB, in addition to DNAext and FM.
% 241010: added beta implementation
% 250422: adding downsampling for fleezers data --for output data only.

% Sample usage:
% [DNAext, FM, fake_offset, FA, FB]=return_FEC('220701', 132,133, 141); 


function [DNAext_clean, FM_clean, fake_offset, FA_clean, FB_clean]= return_FEC(Date,calfilenumber, offsetfilenumber, fextnumber, ...
    manfitrange, usetrapF, beaddiams, method, window)
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
    beaddiams= calparams(1:2);
    method = 'moving'; window = [];
end

if nargin <6
    usetrapF = 'AB'; % by default, use forces from both traps
    beaddiams= calparams(1:2);
    method = 'moving'; window = [];
end

if nargin <7
    beaddiams= calparams(1:2);
    method = 'moving'; window = [];
end

if nargin <8
   method = 'moving'; window = [];
end

tempcalparams = calparams;
tempcalparams(1) = beaddiams(1);
tempcalparams(2) = beaddiams(2);

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
downsamp_fl_data = 1; % downsample fleezers data to ~100 Hz? 0 -- no; 1 -- yes.

param.Date = Date;
param.calfilenumber = calfilenumber;
param.fextnumber = fextnumber;

param.offset_file = 1; % 1 - use offset file, 0 - subtract only QPD offset
if param.offset_file ==0
    warning('Only QPD offset is being subtracted -- offset file not used.');
end
param.file = [param.Date fbslash '_' num2str(param.fextnumber,'%03d')];

startpath = [datadirectories param.Date fbslash];

% calibrate
param.plotCal = 0; % 1 - plot, 0 - do not plot
calRootFile = [param.Date '_' num2str(param.calfilenumber,'%03d') '.dat']; 
calfilename = ['cal' calRootFile(1:(end-4))];
calData = Calibrate(param.Date,calRootFile,tempcalparams,param.plotCal, 0, Savefigdir);


% get raw data
fextRootFile = [param.Date '_' num2str(param.fextnumber,'%03d') '.dat']; 
fextfilename = ['data_' fextRootFile(1:(end-4))];
fextData = ReadJointFile(startpath,fextRootFile);
param.instrument = fextData.oldtrap ; % 1 - old-trap, 0 - fleezers
if fextData.scandir == 0
param.direction = 'Y'; % Y or X
elseif fextData.scandir == 1
param.direction = 'X'; % Y or X
end

if param.instrument
    % get fixedxy results from bead grid analysis -- for old trap
    beta = getBetaFactors(startpath, param.Date); disp('Using beta from file.')
else
    beta = eye(2); % use dummy identity beta for fleezers
end


% get offset data
% offsetfilenumber = param.calfilenumber + 1;
offsetRootFile = [param.Date '_' num2str(offsetfilenumber,'%03d') '.dat']; 
offsetfilename = ['offset' offsetRootFile(1:(end-4))];
offsetData = ReadJointFile(startpath, offsetRootFile);  

if downsamp_fl_data && ~param.instrument  % downsample fleezers data (with averaging)
    fextData=fl_filter_decimate(fextData);
    offsetData=fl_filter_decimate(offsetData);
    disp('Downsampling fleezers data for FEC.');
end

fextData = trap_data_nmConversion(fextData, param, beta); % convert from V or MHz to nm
offsetData = trap_data_nmConversion(offsetData, param, beta); % convert from V or MHz to nm
% subtract offset

if strcmp(usetrapF, 'AB')
    trapData = force_offset_subtraction(offsetData, fextData, calData, param, method, window);  % using XY offset subtraction code
else
    error('Single-trap force offset subtraction not currently enabled.')
    %trapData = force_offset_subtraction_1trap(offsetData, fextData, calData, param, usetrapF);
    %disp(['Manual trap separation correction: ' num2str(trap_sep_corr) 'nm'])
end


%% load construct parameters

ConstructConst = ConstructConstants(construct);

if ~isempty(manfitrange)
    ConstructConst.fitrange = manfitrange;
end

%% calculate but do not plot the WLC model
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
else
    WLC_ext = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp;
end


%% calculate but do not plot the data
% Filter TrapSep
%TrapSepFilt  = movmean(trapData.TrapSep,30);
TrapSepFilt  = movmean(trapData.TrapSep,31); % window should be odd, ideally (AT)
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

%Added 191218, AVT
if isempty(numMAX)
    [numMAX, indMAX]=max(TrapSepFiltClean);
end

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
% Extension limit added 200820, to remove low-extension sticking from fits (AVT)
DNAext_offset = mean(DNAext(FM > ConstructConst.fitrange(1) & FM < ConstructConst.fitrange(2) & DNAext>=150));
WLC_offset = mean(WLC_ext(F > ConstructConst.fitrange(1) & F < ConstructConst.fitrange(2)));
fake_offset_0 = DNAext_offset - WLC_offset;

% SY modified code. Will comment above until SM meeting 012119 
fake_offset = getFakeoffset(DNAext, FM, ConstructConst.fitrange, fake_offset_0, construct) ;

% substract the fake offset
DNAext = DNAext-fake_offset;

%Added on 191112 (AT)
DNAext_clean = [];
FM_clean = [];
FA_clean = [];
FB_clean = [];

for ii = 1:1:rasterNUM
    indxS = indMIN(ii);
    indxE = indMAX(ii);
    indxEE = indMIN(ii+1);
    
    DNAext_clean = [DNAext_clean DNAext(indxS:indxE) DNAext((indxE+1):indxEE)];
    FM_clean = [FM_clean FM(indxS:indxE) FM((indxE+1):indxEE)];
    FA_clean = [FA_clean FA(indxS:indxE) FA((indxE+1):indxEE)];
    FB_clean = [FB_clean FB(indxS:indxE) FB((indxE+1):indxEE)];
    
end

if round(rasterNUM)-rasterNUM == 0.5
    if isempty(ii)
        ii=1;
    end
     indxS = indMIN(ii);
    indxE = indMAX(ii);
    DNAext_clean = [DNAext_clean DNAext(indxS:indxE)];
    FM_clean = [FM_clean FM(indxS:indxE)];
    FA_clean = [FA_clean FA(indxS:indxE)];
    FB_clean = [FB_clean FB(indxS:indxE)];
end


end