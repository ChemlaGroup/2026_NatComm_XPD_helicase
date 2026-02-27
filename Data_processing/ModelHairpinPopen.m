function [P_open,Pn, Hairpin_Position] = ModelHairpinPopen(Seq_num,loop,F,n)

% F is the applied force
% loop is the basepairs at the tight turn on the hairpin. 
% if Seq_num is an integer (#1-6) this code will use existing sequences we have oligos for (see below for seqs).
% if Seq_num is a vector, this code will treat it as the hairpin sequence.
% added 3/7/17: n is # of base pairs 

if nargin <2
    loop = 'TTTT'; F = 12; n=1;
elseif nargin <4
    n=1;
end
if nargin == 0;
Seq_num = 5;
end;

if length(Seq_num) ==1
    switch Seq_num
        case 1 % Zhi HP1
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
        case 8 % HP 8
            seq = 'GA GTC AGT CTC AGT CAC TCA TGT CAG TCA CAG TCA GAG TCA TGT CTG AGT CTT GAT GAT GTC ACT GAC TGA GAC TCT GAC TCT GGC';
        case 9 %HP9
            seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGT CTT GAT GAT GTC ACT GAC TGA GAC TCT GAC TCA CTG AGT CAT GTC TGA GAG TCG';
        case 10 % HP 7.3, T47A, aka HP10
            seq = 'AG TCT CAG TCA CTC ATG TCA GTC ACA GTC AGA GTC ATG TCT GAG TCA TGA TGA TGT CAC TGA CTG AGA CTC TGA CTC ACT GAG TCG';
        case 11 % HP 11 -- 1st, 3rd TGT motifs modified
            seq = 'AGTCTCAGTCACTCATGACAGTCACAGTCAGAGTCATGTCTGAGTCTTGATGACAGTCATGACTGAGACTCTGACTCACTGAGTCG';
        case 12 % HP 12
            seq = 'CAGTCTGAGTCAGTCAGAGTCTGTCAGAGTCAGTCAGAGTCAGTCAGAGTCAGTCAGAGTCAGTCAGAGTCACTGACACTGAGCCG';
        case 13 % HP 13
            seq =  'GTCACTCAGAGACTCACACTCTTGTCATGAGAGTGTTGAGAGAGAGTGAGTGTCTGTGAGACTTGATGTGAGTGTCAGACAGCCGC';
    end
else
    seq = Seq_num;
end
    idx = findstr(seq,' '); seq = seq(setxor(1:length(seq),idx)); %remove spaces

l = length(seq)+length(loop)/2;
% Constants
kT = 1.381e-2*298; %[pN nm]

Pss = 1.0;
Sss = 1000;
hss = 0.6;
lss = length(seq)+length(loop);

if F == 0
    Gss = 0;
else    
    Gss = quad(@(x)XWLCContour(x,Pss,Sss),0,F);
end;

for i = 1:lss+1
    Gss_mat(i) = Gss*hss*2*(i-1);
end;

dGDNA = SeqFreeEnergy2(seq,loop,'FR',0.020 + 0.1 + 8.0*sqrt(0.003))*kT;
GDNA_mat = cumsum([0 dGDNA]); 

Gtot = -GDNA_mat - Gss_mat;
Z = sum(exp(-Gtot/kT));
Pn = exp(-Gtot/kT)./Z; %Prob that n base pairs are open at force F

for i = 1:lss+1
    Hairpin_Position(i) = i-1;
%     P_open(i) = sum(Pn(i+1:l))./sum(Pn(i:l)); %Prob that 1 or more bp are open at position i
    P_open(i) = sum(Pn(i+n:l))./sum(Pn(i:l)); %Prob that n or more bp are open at position i
end;




