function [F_av,xtot,Gtot] = ModelHairpinUnwind(Seq_num,loop,ssbs)

% Created July 11, 2017 - BKS
% Calculates free energy and theoretical model of hairpin pulling curve
% Calls XWLCCountor and SeqFreeEnergy2. Adapated from ModelHairpinUnwindKW
% to include an ssDNA loading site in pulling curve calculation
% 
% Inputs: Seq_num  -  either a number corresponding to a saved sequence
% (below in case structure) or the full 5'->3' sequence of your choice.
% loop  -  hairpin loop sequence
% ssbs  -  length of single stranded binding/loading site

% Updated 210408 by AVT for use with WLC global parameters as declared in
% config file

global WLC_param


% Input hairpin sequence
if nargin <2
    loop = 'TTTT';
    ssbs = 10; %10nt binding site
end
if nargin == 0
Seq_num = 1;
end

if length(Seq_num) ==1
    switch Seq_num
        case 1 % Zhi's HP1
            seq = 'GGC TGA TAG CTG AGC GGT CGG TAT TTC AAA AGT CAA CGT ACT GAT CAC GCT GGA TCC TAG AGT CAA CGT ACT GAT CAC GCT GGA TCC TA';
        case 2
            seq = 'AAT CAG CGA TCA GAT AAC TAA CGC CCT GGG GAC TGG TAC GTC AGC TGT ATC AAG CTT CGA GAC TGG TAC GTC AGC TGT ATC AAG CTT CG';
        case 3
            seq = 'AAT CAT CGA TAA TAT AAC TAA TGC ATC GAA AAT CAG TGA AAA TCA GCT ACA ACG CCC TGG GGA CTG GTA CGT CAG CTG TAT CAA GCT TA';
        case 4 % BKS 1.0
            seq = 'TCT CTC AGT CTG AGA TGT CAG TCT CAG TCA GAG TCT TGT CTG TGT CTT GTT GAT GTC ACT GAC TGA GAC TCT GAC TGT CTG GGT CGC GC';
        case 5 % BKS 1.1
            seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGA GTC ATG TCT GAG TCT TGA TGA TGT CAC TGA CTG AGA CTC TGA CTC ACT GAG TCG AGC';
        case 7 % HP 7.3
            seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGA GTC ATG TCT GAG TCT TGA TGA TGT CAC TGA CTG AGA CTC TGA CTC ACT GAG TCG';
        case 8 % HP 8 -- based on HP 7.3, shifted by 5bp
            seq = 'GA GTC AGT CTC AGT CAC TCA TGT CAG TCA CAG TCA GAG TCA TGT CTG AGT CTT GAT GAT GTC ACT GAC TGA GAC TCT GAC TCT GGC';
        case 9 % HP 9 -- based on HP 7.3, interior portions swapped
            seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGT CTT GAT GAT GTC ACT GAC TGA GAC TCT GAC TCA CTG AGT CAT GTC TGA GAG TCG';
        case 10 % HP 10 (aka 7_1T, or 7_T47A): one basepair different from HP 7.3
             seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGA GTC ATG TCT GAG TCA TGA TGA TGT CAC TGA CTG AGA CTC TGA CTC ACT GAG TCG';
        case 11 % HP 11 -- 1st, 3rd TGT motifs modified
            seq = 'AGTCTCAGTCACTCATGACAGTCACAGTCAGAGTCATGTCTGAGTCTTGATGACAGTCATGACTGAGACTCTGACTCACTGAGTCG';
        case 12 % HP 12 
            seq = 'CAGTCTGAGTCAGTCAGAGTCTGTCAGAGTCAGTCAGAGTCAGTCAGAGTCAGTCAGAGTCAGTCAGAGTCACTGACACTGAGCCG';
        case 13 % HP 13
            seq =  'GTCACTCAGAGACTCACACTCTTGTCATGAGAGTGTTGAGAGAGAGTGAGTGTCTGTGAGACTTGATGTGAGTGTCAGACAGCCGC';
        case 18 % HP 18; mismatch control
             seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGT CTC GAT GAT GTC ACT GAC TGA GAC TCT GAC TCA CTG AGT CAT GTC TGA GAG TCG';
        case 19 % HP 19; base construct for 3' modifications; derived from portions of HP 9
            seq = 'TCAGTGAGTCAGAGTCTCAGTCAGTGACATCATCAAGACTGACTGTGACTGACATGAGTGACTGAGACTGTCATGTCTGAGAGTCG';
    end
else
    seq = Seq_num;
end

idx = findstr(seq,' '); seq = seq(setxor(1:length(seq),idx)); %remove spaces

% Constants
kT = 1.381e-2*298; %[pN nm]

% Params
ktrap = 0.3; %[pN/nm]
lds =  3.265e3; %3.05e3; %2e3; %20; %DNA handle total length [bp] 
lss = length(seq)+length(loop)/2; %hairpin length [bp]

Pss =WLC_param.Pss;  %ssDNA persistence length [nm]
Sss =WLC_param.Sss; %ssDNA stretch modulus [pN]
hds = WLC_param.hds;
hss = WLC_param.hss;
Pds = WLC_param.Pds;
Sds = WLC_param.Sds;

Fmax = 20;
%Fmax = 25;
xssmax = (2*lss+ssbs)*hss*XWLCContour(Fmax,Pss,Sss);
xdsmax = lds*hds*XWLCContour(Fmax,Pds,Sds);
xbmax = Fmax/ktrap;
xtotmax = xbmax + xdsmax + xssmax;  

dx = 1; %5 % nm    10^(round(log10(xtotmax))-2);
xtot = 800:dx:xtotmax;

clear F
F = nan(length(xtot),lss+1);
for i = 1:lss+1
    for j = 1:length(xtot)
        F(j,i) = fzero(@(y) (2*(i-1)+ssbs)*hss*XWLCContour(y,Pss,Sss)+lds*hds*XWLCContour(y,Pds,Sds) + y/ktrap - xtot(j), 1);
    end
end

% Unwinding free energy
% dG = SeqFreeEnergy3(seq,loop,'FR', 0.02+0.05, 0.003)*kT; %[pN nm]  % Kevin's salt correction
% dG = SeqFreeEnergy2(seq,loop,'FR',0.020 + 0.1 + 3.3*sqrt(0.003))*kT; % best fit for HP Seq 1 (I think?)
dG = SeqFreeEnergy2(seq,loop,'FR',0.020 + 0.1 + 8.0*sqrt(0.003))*kT; % used in Zhi's paper, best fit for all those seqs
%AT: vectors in Gtot calculation end up different lengths because dG too
%long. Taking off number of terms equivalent to half the loop at the end;
%they're set to 0 anyways.
if length(loop)>0
    dG = dG(1:end-length(loop)/2);
end
dGDNA = meshgrid(cumsum([0 dG]),ones(1,length(xtot)));

% Trap
xb = F/ktrap;
Eb = 0.5*ktrap*xb.^2; %Trap energy

% ssDNA & dsDNA
dGds = [];
dGss = [];
for i = 1:lss+1
    for j = 1:length(xtot)
        dGss(j,i) = (2*(i-1)+ssbs)*hss*(F(j,i)*XWLCContour(F(j,i),Pss,Sss)-integral(@(f)XWLCContour(f,Pss,Sss),0,F(j,i))); %
        dGds(j,i) = lds*hds*(F(j,i)*XWLCContour(F(j,i),Pds,Sds)-integral(@(f)XWLCContour(f,Pds,Sds),0,F(j,i)));
    end
end

Gtot = Eb + dGds + dGss - dGDNA;
F_av = sum(F.*exp(-Gtot/kT),2)./sum(exp(-Gtot/kT),2); F_av = F_av(:)'; %Sum over all intermediate unwound states i
xtot = xtot-F_av/ktrap;

 figure;
plot(xtot,F_av); 