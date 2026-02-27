Set of MATLAB codes necessary for processing and analyzing data measuring the base pair scale dynamics of XPD helicase from Ferroplasma acidarmanus on DNA hairpins.

The code is fully operational on MATLAB R2021a on Ubuntu Linux (24.04.04 LTS and earlier). The code is expected to be compatible with later versions of MATLAB, and is compatible in full or in part with some earlier versions. The code will also run on other operating systems with a suitable version of MATLAB, after adjustment of file path format via the 'fbslash' global variable and the paths listed in the .config file. No non-standard hardware is required, and no installation is necessary beyond installing MATLAB. 

Code in the 'Calibration' and 'Data_processing' folders was used for analyzing and processing raw data files: calibration and offset files, force-extension curves, and constant-force traces. (Most of these codes are shared; written, modified and updated by Chemla lab members over many years.) Code in the 'Trace_based_analysis' folder then is used specifically for analyzing the processed constant-force data to obtain quantities such as mean residence times as a function of DNA fork position. 

The file "sample_traces_abasic_5set.mat" contains several processed constant-force traces of XPD activity to demonstrate the use of the analysis code in 'Trace_based_analysis', as described below. Plots generated from the sample dataset by some of these codes are provided in 'Trace_sample_analysis_output'.

-----------------------------------------------------------------
Instructions for use:

The 'config.m' file should be executed prior to running other code in order to establish working directories, global variables used throughout the code, display preferences, etc. The config file also lists the correspondence between each construct and its numerical index, which is its identification throughout the code. (See definition of 'HPnames' towards the end of the file.) The sample traces are taken from the 5' set abasic dataset, which has numerical construct index 14. 

The sample data can be loaded by:
data = LoadData(0, 14);

'data' is a cell array such that data{x}{y} corresponds to trace number 'y' for the construct indexed by 'x'. Consequently, data{14}{1} refers to the first trace in the sample dataset. Each element, data{x}{y}, is a structure array containing various kinds of data associated with the trace. 

Plotting the DNA fork position as a function of time for a single trace (after XPD helicase activity begins): 
plotTraces(data{14}) 

Plotting individual bursts instead of whole traces: 
plotBursts(data{14})

To calculate mean residence time versus DNA fork position during unwinding:
[~] = HPlifetimes_errors(0, 1, 1, 25, 'Mean', data, [14]);

This code also returns the number of bursts contributing at each position as well as cumulative residence times over an interval of interest.
The second input determines which part(s) of bursts are included in the calculation. A value of 1 toggles the analysis to unwinding only; -1 to rezipping only; and 0 includes the entire burst. 

As described in the manuscript, we restrict the rezipping analysis to bursts displaying gradual rezipping with no more than short slips, indicating XPD translocating down the 3' side of the hairpin. For accuracy, therefore, the rezipping analysis should not be done on all bursts of a dataset. The optional input variables 'burst_list' and 'timed_burst_list' restrict the analysis to full and/or partial rezipping segments, respectively, of the burst subsets (manually) selected to these lists.

If a reference data set for an unmodified construct with the same sequence is available, the excess dwell time on a modified construct can be calculated by the function 'HPlifetimes_errors_differences'.

The probability of traversal past a threshold (default: 30 bp) during unwinding can be found for a single dataset by:
[survival_frac, times, nbursts] = BurstSurvivalFrac(data{14}, 1, 30,[]);

or for multiple datasets in tandem by:
BurstSurvivalFrac_plotDatasets(data, [clist]), where clist is a list of numerical indeces of datasets of interest.

BurstSurvivalFracRezipping_plotDatasets returns traversal probabilities for rezipping. As described in the manuscript, the analysis was carried out on a subset of bursts displaying gradual rezipping, and the input burst_list is constructed from the (manual) selection of this subset.  

Bootstrapping can be done on the residence time distributions to identify outlier bursts:
HPlifetimes_bootstrap(0, 1, 1, 25, 'Mean', data, [14])

The plots of distributions of the means at different hairpin positions can then be assessed for multimodality, indicating the presence of an outlier which skews the distribution. This bootstrapping procedure takes only several seconds to execute on an ordinary computer. All other codes described above run on a similar or faster timescale.    

For outputs of some of the above commands for the sample dataset, please see 'Trace_sample_analysis_output'.
