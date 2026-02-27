% Author: Monika A. Makurath
% Date created: 06/12/2018
% This is a configuration file for MATLAB analysis on the old trap.

% Modified by Alice Troitskaia for analysis-specific configuration.

global datadirectories
global trap_computer
global analysis
global calparams
global fbslash
global WLC_param
global construct
global HPcolours
global HPnames


trap_computer = 0;
fbslash = '/'; % could also do filesep

main_dir = pwd;

analysis = main_dir; % modify as needed; assuming main directory containing code subfolders, added to path below.

           %beadA  beadB    fsamp      fxyhi    fsumhi  flo    avwin
calparams = [845    818     125000      62000     8000   600     200]; disp('OT configuration!')
%calparams = [818    845     125000      33000     8000   600     200]; disp('Fleezers configuration!')
%calparams = [818    845     125000      18000     8000   600     200]; disp('Fleezers configuration! (Spikes)')
% Using higher flo for data with bad calibrations (noisy at low frequencies). With good calibrations,
% can use 600.

% Using same values as Barbara Stekas:
WLC_param.hds = 0.34; % [nm/bp] interphosphate distance 
WLC_param.hss = 0.6; % nm/nt
WLC_param.Pds = 50; % nm
WLC_param.Pss = 1.07; % nm -- changed 200204 from 1.0 nm, to value Steve had arrived at a while ago
WLC_param.Sss = 1000; % pN
WLC_param.Sds = 1000; % pN
WLC_param.kT = 4.14;% Zhi Qi's values (dsDNA)

construct = 'Hairpin-LambdaRH-dT-10-86bp';


datadirectories = [main_dir '/']; % modify as needed

% Generate hairpin unwinding model if wanted later for plotting with FECs
%[HPF_av,HPxtot,hpGtot] = ModelHairpinUnwind(9,'TTTT',10);hpmodel = [HPxtot; HPF_av];

figdir = [main_dir '/']; % modify as needed


% For sequence analysis:

hp_sequences{9} =  'AGTCTCAGTCACTCATGTCAGTCACAGTCAGTCTTGATGATGTCACTGACTGAGACTCTGACTCACTGAGTCATGTCTGAGAGTCG';
hp_sequences{19} = 'TCAGTGAGTCAGAGTCTCAGTCAGTGACATCATCAAGACTGACTGTGACTGACATGAGTGACTGAGACTGTCATGTCTGAGAGTCG';

%Color scheme

HPcolours{9}  = [0.502   0.502      0.502];
HPcolours{14} = [0.8706  0.4902     0];
HPcolours{24} = [0       0.7490     0.7490]; % off-strand mismatch
HPcolours{16} = [1       0          0];
HPcolours{17} = [0       0.4980     0];
HPcolours{15} = [0       0.7490     0.7490];%[0 0 0]; % for plotting if necessary


HPcolours{19} = [0.502   0.502      0.502];
HPcolours{22} = [0.8706  0.4902     0];
HPcolours{27} = [0       0.7490     0.7490]; % off-strand mismatch
HPcolours{20} = [1       0          0];
HPcolours{21} = [0       0.4980     0];
HPcolours{23} = [0       0.7490     0.7490];%[0 0 0]; % for plotting if necessary


HPcolours{1}=[0 0 0]; % dummy; miscellaneous testing purposes
HPcolours{2}= [1 0.07 0.65]; %dummy; miscellaneous testing purposes
HPcolours{3}=[0 1 0]; % dummy; miscellaneous testing purposes


HPnames{9}= '5'' Uniform';
HPnames{14}= '5'' Abasic';
HPnames{15}= '5'' On-strand Mismatch';
HPnames{16}= '5'' CPD';
HPnames{17}= '5'' Fluorescein';
HPnames{19}= '3'' Uniform';
HPnames{20}= '3'' CPD';
HPnames{21}= '3'' Fluorescein';
HPnames{22}= '3'' Abasic';
HPnames{23}= '3'' On-strand Mismatch'; 
HPnames{24}= '5'' Off-strand Mismatch';
HPnames{27}= '3'' Off-strand Mismatch';

HPnames{1}=  'Dummy 1';   % miscellaneous testing purposes
HPnames{2}=  'Dummy 2';   % miscellaneous testing purposes
HPnames{3}=  'Dummy 3';   % miscellaneous testing purposes


% Set plotting defaults (from Barbara, added 2017/02/27) -- modified after
% matlab version upgrade
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultLineMarkerSize', 6);
set(0,'DefaultUicontrolFontSize', 16);
set(0,'DefaultFigureWindowStyle','docked');
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1) % added 240228

% Improving the new and unpleasant Matlab defaults...
set(groot, 'defaultAxesXColor', [0,0,0], ...
           'defaultAxesYColor', [0,0,0], ...
           'defaultAxesZColor', [0,0,0]);
% Reversing font smoothing of newer versions
set(groot,'defaultAxesFontSmoothing', 'off'); 
set(groot,'defaultTextFontSmoothing', 'off'); 

cd(analysis)


addpath(analysis, [analysis fbslash 'Calibration'],[analysis fbslash 'Data_processing'], [analysis fbslash 'Trace_based_analysis']);