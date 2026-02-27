
% Copied on 220906 from earlier separate project

% Plot FECs, colour-code by bead pair -- modified, standalone version of
% section from FEC_overlay_noofset

% Inputs:
%  flist, comprising FEC metadata (Date; Cal; Offset; FEC;	Same
% molecule as last? -- binary	New bead pair? -- binary)

% Buffer -- for plot-labeling purposes

%offset:  0 for no offset subtraction; 1 for subtraction of offset calculated by
% return_FEC; any other value for direct extension shift by that value

% Output:
% fake_offsets_all -- vector of 'fake' offsets; obtained by fitting each
% FEC to theoretical model within a range set in return_FEC (differs for
% different constructs). 

%230329: modified to load construct using global definition for model
%plotting
% 200323: modified to separate modification of DNA_ext and plotting, to
% enable proper rescaling by values at 10pN.

% [fake_offsets_all, ~] = FEC_overlay_processing(flist,1, 'XPD buffer');

function [fake_offsets_all, beaddiams_all] = FEC_overlay_processing(flist, offset, buffer)

global calparams
global WLC_param
global construct
global HPcolours

if nargin<3
%buffer = 'Imaging buffer'; buffer = 'TMS20';
buffer = 'XPD';
offset = 0;
end

%fitrange = [5 15];
%fitrange = [4 10];
fitrange = [5 10];
%fitrange = [];
scale10pNext=0;
scale10pNF =0;
%colour_by_beadpair = 1; % if 0, colour by molecule
colour_by_beadpair = 0; 
optimize_beaddiams=0;

% Enforce first bead/tether in list being new. Sometimes necessary when
% considering subset of list.
flist(1, 7)= 0; % new tether

if colour_by_beadpair
    flist(1, 8)=1; % new bead
end

if colour_by_beadpair
    %Nbeads = nansum(flist(:, 8));
    Nbeads = nansum(flist(:, 8)~=0);
    Nbeadcolours = distinguishable_colors(Nbeads);
end

legendstr = {};

figure; hold on;

% load construct parameters
ConstructConst = ConstructConstants(construct);
if isempty(fitrange)
    fitrange = ConstructConst.fitrange;
end

F_model = 0.31:0.01:20;
bp = 1*WLC_param.hds*XWLCContour(F_model,WLC_param.Pds,WLC_param.Sds,WLC_param.kT); % nm/bp
nt = 1*WLC_param.hss*XWLCContour(F_model,WLC_param.Pss,WLC_param.Sss,WLC_param.kT); % nm/nt
ext_model_closed = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp; %10*nt + 3265*bp;
ext_model_open = ConstructConst.ssDNA2.*nt + ConstructConst.dsDNA2.*bp; %(10+86*2+4)*nt + 3265*bp;

plot(ext_model_closed, F_model, '--r');
plot(ext_model_open, F_model, '--k');
legendstr = [legendstr 'WLC closed' 'WLC open'];


Nmolecules = size(flist, 1)-nansum(flist(:, 7));

Nmoleculecolours = distinguishable_colors(Nmolecules);

bead_count = 0;
molecule_count = 0;
N = size(flist, 1);

fake_offsets_all = [];
beaddiams_all = [];



for i=1:N
    
    
  if colour_by_beadpair  
    if flist(i, 8) ~=0
        bead_count = bead_count+1;       
    end
  else 
      if flist(i, 7) ==0
        molecule_count = molecule_count+1;
    end
  end
    
    if optimize_beaddiams
        disp('Bead diameter optimization!')
        opt_forcerange = [3 14];
        %opt_forcerange= fitrange;
        % Calculate theory to fit to:
        F_model = 0:0.001:20;
        bp = 1*WLC_param.hds*XWLCContour(F_model,WLC_param.Pds,WLC_param.Sds,WLC_param.kT);
        nt = 1*WLC_param.hss*XWLCContour(F_model,WLC_param.Pss,WLC_param.Sss,WLC_param.kT);
        ext_model = 44*nt + 1710*bp;    % Steve's fork-hairpin construct
        %Both beads:
        %bd = fminsearch(@(beaddiams) HPModelFit(beaddiams, F_model, ext_model, num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), 'B', opt_forcerange, fitrange), [2100 818]);
        % ADig bead only:
           bd = fminsearch(@(Abeaddiam) HPModelFit([Abeaddiam calparams(2)], F_model, ext_model, num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), 'B', opt_forcerange, fitrange), 2100);

        disp(bd);
        beaddiams_all = [beaddiams_all bd];
        % Calculate force, extension, offset for optimized bead diameters
    [DNAext, FM, fake_offset] = return_FEC(num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), fitrange, 'B',[bd calparams(2)]);
    
    
    else
        [DNAext, FM, fake_offset] = return_FEC(num2str(flist(i, 1)), flist(i, 2), flist(i, 3), flist(i, 4), fitrange, 'AB',calparams(1:2));
    end
    
    
    
    if isnan(fake_offset)
        disp(['Bad fake offset: ' num2str(flist(i, 1)) ' ' num2str(flist(i, 4))]);
        
    else
        fake_offsets_all = [fake_offsets_all fake_offset];
    end
    
    if flist(i, 1) == 200201 && flist(i, 4) == 45
        DNAext = DNAext(1:end-15);
        FM = FM(1:end-15);
    end
    
    % Offset corrections. If 1, leave as is (i.e. with fake_offset
    % incorporated by return_FEC code.)
     if offset ==0
        % Add the fake offset back in, to get data unshifted in extension
        DNAext = DNAext+fake_offset;
     elseif offset ~=0 && offset ~=1
         DNAext =DNAext+fake_offset + offset;
    end
    
    %scale by 10 pN extension:
    if scale10pNext
        ind10pN = find(abs(FM-10)<0.5);
        % If FEC does not reach 10 pN, skip it and move to the next one.
        if isempty(ind10pN)
            continue
        end
        mean10pN_F=mean(FM(ind10pN));
        mean10pN_ext = mean(DNAext(ind10pN));
        disp(['Scaling by ' num2str(mean10pN_F) ' pN, ' num2str(mean10pN_ext) ' nm' ])
        DNAext = DNAext/mean10pN_ext;
        if scale10pNF
            FM = FM/mean10pN_F;
        end
    end
    
% Modify extension earlier.    
%     if offset ==0
%         % Add the fake offset back in, to get data unshifted in extension
%         plot(DNAext+fake_offset, FM, 'Color', Nbeadcolours(bead_count, :));
%     elseif offset ==1
%         plot(DNAext, FM, 'Color', Nbeadcolours(bead_count, :));
%     else
%         plot(DNAext+fake_offset + offset, FM, 'Color', Nbeadcolours(bead_count, :));
%     end
    
    % Actual plotting
    if colour_by_beadpair
        Ncolor = Nbeadcolours(bead_count, :);
    else
        Ncolor = Nmoleculecolours(molecule_count, :);
    end
    
    %if specific colour needed
    %Ncolor = HPcolours{27};
    
    %plot(DNAext, FM, 'Color', Nbeadcolours(bead_count, :));
    plot(DNAext, FM, 'Color', Ncolor);
    
    %legendstr{i} = [num2str(flist(i, 1)) '\_' num2str(flist(i, 4))];    
    legendstr = [legendstr [num2str(flist(i, 1)) '\_' num2str(flist(i, 4))]]; 
end

av_fake_offset = sum(fake_offsets_all)/N;
disp(['Average fake offset from fit to WLC: ' num2str(av_fake_offset)]);

legend(legendstr)
xlabel('Extension (nm)');
ylabel('Force (pN)');
%title('Force vs. Extension of ssDNA (1180nt)')
if offset ==0
    %title([{'Force vs. Extension of gxDNA (1180nt); no extension offset'}
      %  {[buffer ', ' num2str(flist(1, 1))]}]);
    figtitle = [{'Force vs. Extension; no extension offset'}
        {[buffer ', ' num2str(flist(1, 1))]}];
elseif offset ==1
    
    %title([{'Force vs. Extension of gxDNA (1180nt); aligned to mean WLC value in range (5, 15) pN'}
       % {[buffer ', ' num2str(flist(1, 1))]}]);
    figtitle= [{['Force vs. Extension; aligned to mean WLC value in range (' num2str(fitrange(1))  ', ' num2str(fitrange(2)) ') pN']}
        {[buffer ', ' num2str(flist(1, 1))]}];
    
else
     %title([{['Force vs. Extension of gxDNA (1180nt); extension offset ' num2str(offset) ' nm']}
        %{[buffer ', ' num2str(flist(1, 1))]}]);
    figtitle= [{['Force vs. Extension; extension offset ' num2str(offset) ' nm']}
        {[buffer ', ' num2str(flist(1, 1))]}];
    
end

if scale10pNext
    figtitle{end+1} = 'Scaled by 10 pN values';
    xlabel('L/L_{10 pN}');
    if scale10pNF
        ylabel('F/F_{10 pN}');
    end
end
 
 figtitle{end+1} = ['Average fake offset of fits to WLC: ' num2str(av_fake_offset, '%.1f') ' nm'];
 
 title(figtitle)
 legend('off')
 % To reverse axes and view as log-log plot:
 %view(-90, 90)
 %set(gca, 'ydir', 'reverse');
 %ylim([0.01 35])
 %set(gca, 'YScale', 'log')
 %set(gca, 'XScale', 'log')
 
 
 end