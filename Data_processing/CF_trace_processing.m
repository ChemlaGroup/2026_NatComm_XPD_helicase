
% AVT
% Code for processing constant-force data traces, with XPD activity.

% Based on Barbara's Fluorescence_HP_Analysis and Monika's
% plot_force_extension (from the old trap).
% Uses explicit variables as defined in Monika's code, and as far as
% possible, same variable names.

% Currently not including fluorescence processing, which Barbara's code
% did.

% 230912: modifying to allow data reprocessing: instead of dataSet
% structure, load in an already-processed dataset.

function output=CF_trace_processing()


global datadirectories
global analysis
global calparams
global fbslash


saving = 1;  %Saving figures
Savefigdir = [datadirectories 'Processed/'];
if ~exist(Savefigdir, 'dir')
    mkdir(Savefigdir);
end

%Format {'date' calNumber tetherNumber fileNumber beadsizes selectedStart5Time selectedEndTime usetrapF}

dataSet =  {...

% {'241221' 1 4 5 'AB' [118 121] 3}...
% {'241221' 11 13 14 'AB' [112.5 112.9] -1}...
% {'241221' 40 42 43 'AB' [212 217] 3}...
% {'241221' 60 67 68 'AB' [150 160] -1}...
% {'241221' 60 70 71 'AB' [150 170] -1}...
% {'241221' 60 72 73 'AB' [115 125] -1}...
% {'241221' 60 74 75 'AB' [130 150] -1}...
% {'241221' 60 76 77 'AB' [149 150.5] -1}...
% {'241221' 60 78 79 'AB' [0 10] -1}...
% {'241222' 21 25 26 'AB' [119 125] 3}...

};

reprocess = 1;

if reprocess ==1
    if ~isempty(dataSet)
        error('Dataset entered but reprocessing mode requested. Verify inputs.')
    end
    %load in an already-processed dataset
    %[dataSet, ~] =  DatasetLoading_Baselineselect(9, 0);

    %load('datadirectories/data241007.mat');
    %dataSet= data241007;
else
    if isempty(dataSet)
        error('No dataset entered. Verify inputs.')
    end
end

downsamp_fl_data = 1; % downsample fleezers data to ~100 Hz? 0 -- no; 1 -- yes.
if downsamp_fl_data
    disp('Fleezers data will be downsampled to ~100 Hz.')
end

output = cell(1,numel(dataSet));

same = 0;
for k=1:numel(dataSet)

    %----------------------------------------------------------------------
    % Read in dataSet info
    if reprocess ==1
        Date = dataSet{k}.date;
        calfilenumber= dataSet{k}.refnums(1);
        fextnumber = dataSet{k}.refnums(3);
        CFnumber= dataSet{k}.refnums(4);
        beadsize = dataSet{k}.cal.beadsize;
        usetrapF = dataSet{k}.usetrapF;
    else
        Date    = dataSet{k}{1};
        calfilenumber    = dataSet{k}{2};
        fextnumber    = dataSet{k}{3};
        CFnumber    = dataSet{k}{4};
        %beadsize = dataSet{k}{5};
        beadsize = calparams(1:2); % [bead A, bead B]
        usetrapF = dataSet{k}{5};
    end
    
    startpath   = [datadirectories Date fbslash];
    CFRootfile    = [Date '_' num2str(CFnumber,'%03d') '.dat'];
    calRootFile     = [Date '_' num2str(calfilenumber,'%03d') '.dat'];
    offsetRootFile     = [Date '_' num2str(calfilenumber+1,'%03d') '.dat'];
    
    output{k}.date = Date;
    output{k}.refnums = [calfilenumber, calfilenumber+1, fextnumber, CFnumber, 1];
    info = dir([startpath CFRootfile]);
    output{k}.usetrapF = usetrapF;

    if reprocess ==1
        output{k}.datevec = dataSet{k}.datevec;
        output{k}.baselines = dataSet{k}.baselines;
        if isfield(dataSet{k}, 'tether_rating')
            output{k}.tether_rating = dataSet{k}.tether_rating;
        end
    else
         output{k}.datevec = datevec(info.datenum);
        if ~isempty(dataSet{k}{6})
            output{k}.baselines = dataSet{k}{6};
        end
        if ~isempty(dataSet{k}{7})
            output{k}.tether_rating = dataSet{k}{7};
        end
    end    
    
    % Get experimental conditions
    if reprocess == 1
        output{k}.conc = dataSet{k}.conc;
        output{k}.RPA2label = dataSet{k}.RPA2label;
        output{k}.buffer = dataSet{k}.buffer;
        output{k}.construct = dataSet{k}.construct;
        output{k}.scav = dataSet{k}.scav;
        output{k}.BSA = dataSet{k}.BSA;
        output{k}.ATPtop = dataSet{k}.ATPtop;
        %output{k}.hp_params = dataSet{k}.hp_params;
        
    else
        if ~same
            prompts = {'RPA2 conc (nM)','XPD conc (nM)','ATP conc (uM)','ATPgammaS conc (uM)','RPA2 label','buffer','construct','ox.scav.','BSA','ATP on top?','Same for rest of data sets?'};
            defaults = {'0','20','500','50','NA','TMS20','10T_HP27','poxy','15','0','1'};
            
            conditions = inputdlg(prompts,['Enter experimental conditions for ' CFRootfile],1,defaults);
            same = str2double(conditions{end});
        end
        
        output{k}.conc = { 'RPA2 nM' 'XPD nM' 'ATP uM' 'ATPgammaS uM' ; str2double(conditions{1}) str2double(conditions{2}) str2double(conditions{3}) str2double(conditions{4})};
        output{k}.RPA2label = conditions{5};
        output{k}.buffer = conditions{6};
        output{k}.construct = conditions{7};
        output{k}.scav = conditions{8};
        output{k}.BSA = conditions{9};
        output{k}.ATPtop = conditions{10};
        
    end
        
        switch output{k}.construct
            case {'10T_HP1','10T_BKS1.1','10T_BKS1.0'}
                hp_params(1) = 3050; %length of ds handles
                hp_params(2) = 10; %length of ss loading site
                hp_params(3) = 89; %length of hairpin
                hp_params(4) = 4; %length of T loop in hairpin
            case 'Cy3_HP'
                hp_params(1) = 3050; %length of ds handles
                hp_params(2) = 10; %length of ss loading site
                hp_params(3) = 50; %length of hairpin
                hp_params(4) = 4; %length of T loop in hairpin
            case {'10T_HP7.3', '10T_HP8', '10T_HP9', '10T_HP10', '10T_HP11', '10T_HP12', '10T_HP13', ...
                    '10T_HP14','10T_HP15','10T_HP16', '10T_HP17', '10T_HP18', '10T_HP19', '10T_HP20',...
                    '10T_HP21', '10T_HP22', '10T_HP23', '10T_HP24', '10T_HP27'}
                hp_params(1) = 3265; %length of ds handles
                hp_params(2) = 10; %length of ss loading site
                hp_params(3) = 86; %length of hairpin
                hp_params(4) = 4; %length of T loop in hairpin
        end
        
        output{k}.hp_params = hp_params;
        
  
    %----------------------------------------------------------------------
    %Calibration; calculate, plot FEC
    
    calData = Calibrate(Date,calRootFile,calparams,0,0);
    
    [DNAext_fext, FM_fext, fake_offset, FA_fext, FB_fext]= return_FEC(Date,calfilenumber, (calfilenumber+1), fextnumber, [], usetrapF);
    output{k}.tether.ext = DNAext_fext; % already incorporates fakeoffset
    output{k}.tether.force = FM_fext;
    output{k}.tether.fakeoffset = fake_offset;

    % Basic plot of FEC with model:
    plot_precalculated_FEC(DNAext_fext, FM_fext, FA_fext, FB_fext); 
    pull = gcf; % FEC figure
    title(['FEC for ' Date '\_' num2str(fextnumber, '%03d') '; fakeoffset ' num2str(fake_offset) ' nm']);
   
    
    %----------------------------------------------------------------------
    % Read data files; save raw data and offset
    CFData = ReadJointFile(startpath,CFRootfile);  
    offsetData = ReadJointFile(startpath, offsetRootFile);
    
    if downsamp_fl_data && ~CFData.oldtrap  % downsample fleezers data (with averaging)
        CFData=fl_filter_decimate(CFData);
        offsetData=fl_filter_decimate(offsetData);
        disp('Downsampling fleezers data for CF.');
        output{k}.rawdata.downsampfactor = CFData.downsampfactor;
    end
    
    
    output{k}.cal.beadsize = beadsize;
    output{k}.cal.cal = calData;
    output{k}.rawdata.data = CFData;
    output{k}.rawdata.offset = offsetData;
    output{k}.oldtrap = CFData.oldtrap;
    % adding direction field, following Monika's plot_force_extension code.
    if abs(CFData.t1x-CFData.t2x) < 0.1
        output{k}.direction = 'Y'; % Y or X
    elseif abs(CFData.t1y-CFData.t2y) < 0.1
        output{k}.direction = 'X'; % Y or X
    else
       disp('Direction not clear, skipping.')
    end
    
    % Read times of stage moves:
    [output{k}.stgTime, output{k}.stgPos] = StageMoveTime(Date, [Date '_' num2str(CFnumber,'%03d') '_stagetime.txt']);
    
    %----------------------------------------------------------------------
    % Calculate extension and force
    param.instrument = CFData.oldtrap;  % 1 - Old Trap, 0 - fleezers
    param.offset_file = 1; % 1 - use offset file, 0 - subtract only QPD offset
    param.Date = Date;
    
    if param.instrument
        % get fixedxy results from bead grid analysis -- for old trap
        beta = getBetaFactors(startpath, param.Date); disp('Using beta from file.')
    else
        beta = eye(2); % use dummy identity beta for fleezers
    end
    
    CFData=trap_data_nmConversion(CFData, param, beta); % convert from V or MHz to nm
    offsetData = trap_data_nmConversion(offsetData, param, beta); % convert from V or MHz to nm

    % use 2D offset subtraction:
    CFtrapdata = force_offset_subtraction(offsetData, CFData, calData, param); % using XY offset
   
    % % Taking both x and y into account
    FA = sqrt(CFtrapdata.force_AX.^2 + CFtrapdata.force_AY.^2);
    FB = sqrt(CFtrapdata.force_BX.^2 + CFtrapdata.force_BY.^2);
    FM = (FA + FB)./2;
    % Baseline subtraction is done via manually-selected window later, but
    % subtracting fake_offset here should give about correct extension and
    % keep bp conversion more accurate:
    CFtrapdata.DNAext=CFtrapdata.DNAext-fake_offset; 
    
    %----------------------------------------------------------------------
    %Select Average Force region
    
    
    % select the force opposite to pulling orientation:
    if strcmp(output{k}.direction,'X')
        CFtrapdata.scandir = 1; % not actually field read in for this type of data
        Apos_orthogonal = CFtrapdata.A_Ynm;
        Bpos_orthogonal = CFtrapdata.B_Ynm;
    elseif strcmp(output{k}.direction,'Y')
        CFtrapdata.scandir = 0; % not actually field read in for this type of data
        Apos_orthogonal = CFtrapdata.A_Xnm;
        Bpos_orthogonal = CFtrapdata.B_Xnm;
    end
    
    
    if reprocess ==1 && isfield(dataSet{k}, 'force_region')
        output{k}.force_region = dataSet{k}.force_region;
    else
        ax = [0 0 0];
        figure;
        subplot(3,2,[1,2]); plot(CFtrapdata.time,CFtrapdata.DNAext); ylabel('Ext (nm)')
        ax(1) = gca; title('Select Average Force region');
        
        subplot(3,2,[3,4]); plot(CFtrapdata.time,FM); ylabel('Force (pN)')
        ax(2) = gca;
       
        
        subplot(3,2,[5,6]); hold on;
        a = plot(CFtrapdata.time, Apos_orthogonal, 'r'); set(a,'DisplayName', 'A position \_orthogonal');
        b = plot(CFtrapdata.time, Bpos_orthogonal, 'b'); set(b,'DisplayName', 'B position \_orthogonal');
        xlabel('Time (s)')
        ylabel(['Bead Position orthogonal' newline 'to pulling direction (nm)'])
        ax(3) = gca;
        
        linkaxes(ax,'x');
        %pause
        ft = ginput(2);
        output{k}.force_region = ft(1:2);       
        close(gcf);
    end
    
    %disp([output{k}.force_region]);
     fi = CFtrapdata.time<= output{k}.force_region(2) & CFtrapdata.time>= output{k}.force_region(1);
    
    %----------------------------------------------------------------------
    % Build output
    
   
    output{k}.time = CFtrapdata.time;
    output{k}.ext = CFtrapdata.DNAext;
    output{k}.traprate = 1/(CFtrapdata.time(2)-CFtrapdata.time(1));
    output{k}.force = FM;
    output{k}.force_A = FA;
    output{k}.force_B = FB;
    output{k}.avgforce = mean(FM(fi), 'omitnan'); 

    % Overtop of FEC
    figure(pull); plot((output{k}.ext),output{k}.force,'k', 'DisplayName', 'CF data');
    
    %----------------------------------------------------------------------    
    % Convert nm to bp (time-dependent) 
 
    CFtrapdata.bpraw = bpunwound(output{k}.force,output{k}.ext,output{k}.avgforce,output{k}.hp_params);
    output{k}.bpraw = CFtrapdata.bpraw;
    
   
    %----------------------------------------------------------------------     
    % Baseline corrections
    % Select points on bp baseline to fit

    if ~isfield(output{k}, 'tether_rating')
        output{k}.tether_rating = -1;
    end
    
    if isfield(output{k}, 'baselines')
        regions = output{k}.baselines;        
    else
        figure; plot(output{k}.time, output{k}.bpraw );
        title('Zoom in find baseline regions at start and end of trace. Press space when ready to enter time windows.');
        xlabel('Time (sec)'); ylabel('Base pairs opened')
        if max(output{k}.bpraw) > 95
            ylim([-5 95])
        elseif min(output{k}.bpraw) <-10
            ylim([-5 max(output{k}.bpraw)+5])
        end
        
        %pause
        
        prompts = {'Early region start', 'Early region end', 'Later region start', 'Later region end'};
        regions = getnumbers('Enter start and end times of desired region:', prompts, {'0','5','50','60'});
        
        output{k}.baselines = regions;
        close(gcf)
        
    end
    
  % Carry out baseline fit  
    
  if length(regions)>2
      if (regions(1)==0) && (regions(2)==0)
          x1 = 0; x2 = 1; y1 = 0; y2=0;
      else
          fitwin1 = find((output{k}.time >= regions(1) & output{k}.time <=regions(2)));
          fitwin2 = find((output{k}.time >= regions(3) & output{k}.time <=regions(4)));
          
          x1 = output{k}.time(fitwin1);
          x2 = output{k}.time(fitwin2);
          y1 = output{k}.bpraw(fitwin1);
          y2 = output{k}.bpraw(fitwin2);
      end
      
      extdrift = (mean(y2)-mean(y1))*(output{k}.time-mean(x1))/(mean(x2)-mean(x1))+mean(y1);
  else
      y1 = output{k}.bpraw(output{k}.time >= regions(1) & output{k}.time <=regions(2));
      extdrift = mean(y1);
  end
  
  output{k}.bp = output{k}.bpraw - extdrift;
    

    %----------------------------------------------------------------------   
    % If reprocessing, save prior burst selections.
    
    if reprocess ==1 && isfield(dataSet{k}, 'burst2')
        output{k}.burst2 = dataSet{k}.burst2;
    end
    
   
    %---------------------------------------------------------------------- 
    % Plots
    
    sp=3; % could make more subplots if have more fields, e.g. fluorescence
    figure; set(gcf,'Name',[Date '_' num2str(CFnumber,'%03d') ' ' num2str(output{k}.traprate) 'Hz']);
    axeshandles = [];
    
    subplot(sp,1,1);
    plot(output{k}.time, output{k}.bp, 'Color', [0.4 0.4 0.4]);
    hold on;
     plot(FilterAndDecimate(output{k}.time,3), FilterAndDecimate(output{k}.bp,3), 'k','LineWidth',1); 
     % Should this use a different filter? This doesn't treat edges well.
     % But that's negligible here.
     
    ylabel('Basepairs Open');
        ylim([-2,100]);
        axeshandles = [axeshandles gca];
      title([Date '\_' num2str(CFnumber,'%03d') ' ' num2str(output{k}.traprate) 'Hz, force from trap(s) ' usetrapF]);    
      
      %Force
      subplot(sp,1,sp-1); hold on;
      plot(output{k}.time, output{k}.force_A, 'Color', [1 0.5 0.5], 'DisplayName', 'Trap A');
      plot(output{k}.time, output{k}.force_B, 'Color', [0.5 1 0.5], 'DisplayName', 'Trap B');
      plot(output{k}.time, output{k}.force, 'Color', [0.4 0.4 1], 'DisplayName', 'Average (A and B)');
      hold on;
      ylabel('Force (pN)')
      if any(isnan(output{k}.force))
          ylim([0 20]);
      else
          ylim([min(output{k}.force)-1,max(output{k}.force)+1]);
      end
      axeshandles = [axeshandles gca];
      
      %Position in direction orthogonal to tension
      
      % select the force opposite to pulling orientation: 
      % need to re-calculate after corrections
    if CFtrapdata.scandir == 1
        Apos_orthogonal = CFtrapdata.A_Ynm;
        Bpos_orthogonal = CFtrapdata.B_Ynm;
        dispnameA = 'A_Y'; dispnameB = 'B_Y';
        ylabel_text ='Bead Position in Y (nm)';
    elseif CFtrapdata.scandir == 0
        Apos_orthogonal = CFtrapdata.A_Xnm;
        Bpos_orthogonal = CFtrapdata.B_Xnm;
        dispnameA = 'A_X'; dispnameB = 'B_X';
        ylabel_text ='Bead Position in X (nm)';
    end
      
      
        subplot(sp,1,sp); hold on;
        a = plot(output{k}.time, Apos_orthogonal, 'r'); set(a,'DisplayName', dispnameA);
        b = plot(output{k}.time, Bpos_orthogonal, 'b'); set(b,'DisplayName', dispnameB);
        xlabel('Time (s)')
        ylabel(ylabel_text);
        axeshandles = [axeshandles gca];
        
        linkaxes(axeshandles,'x')
      
        
        if saving
            
            var = ['tr' Date '_' num2str(CFnumber,'%03d')]; %make trace variable name
            eval([var '= output{k};'])            
            cd(Savefigdir)
            save([Date '_' num2str(CFnumber,'%03d') ' struct.mat'], var);
            saveas(gcf,[Date '_' num2str(CFnumber,'%03d') '_'  num2str(output{k}.traprate) 'Hz_' output{k}.construct '_' num2str(output{k}.conc{2,1}) 'R_' num2str(output{k}.conc{2,2}) 'X_' num2str(output{k}.conc{2,3}) 'A.fig']);
            saveas(pull,[Date '_' num2str(CFnumber,'%03d') '_FEC.fig']);
            cd(analysis)
        end

        
      
end

end

