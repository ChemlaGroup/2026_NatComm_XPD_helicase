function [ConstructConst] = ConstructConstants(ConstructName)
% 190116 Steve Yeo 
% Output ConstructConst constants like: 
% hairpin : Is the contsruct hairpin or not ( 1: T , 0: F ) 
% dna1label = name of dna1 label 'Hairpin Closed Model';
% ssDNA1 : length of ssDNA1 in closed configuration 
% dsDNA1 : length of dsDNA1 in closed configuration 
% dna2label = name of dna2 label 'Hairpin Closed Model';
% ssDNA1 : length of ssDNA in open configuration 
% dsDNA1 : length of dsDNA in open configuration 
% fitoffset : ??? 
% fitclosed : ??? 
% fitrange : range of force to fit model and FEC 
% fakeoffset : ??? 
% stdE: ??? 
% stdE2: ??? 



switch ConstructName
    case '105dt-SSB-gapped'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 105; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710+1555;
%         ConstructConst.dna2label = 'Hairpin Open Model';
%         ConstructConst.ssDNA2 = 4+20+1719*2; %2*50+4+10; %
%         ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 6]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
    
    
    case '35dt-SSB-gapped'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 35; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710+1555;
%         ConstructConst.dna2label = 'Hairpin Open Model';
%         ConstructConst.ssDNA2 = 4+20+1719*2; %2*50+4+10; %
%         ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
    
    
    case '5Fork-Hairpin-20dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+20; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+20+1719*2; %2*50+4+10; %
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case '5Fork-Hairpin-10dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+10; %4dT spacer and 10dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+10+1719*2; 
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case '5Fork-Hairpin-30dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+30; %4dT spacer and 10dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+30+1719*2; 
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
        
        
        
    case '5Fork-Hairpin-40dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+40; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+40+1719*2; %2*50+4+10; %
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
        
   case '5Fork-Hairpin-50dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+50; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+50+1719*2; %2*50+4+10; %
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case '5Fork-Hairpin-60dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+60; %4dT spacer and 20dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+60+1719*2; %2*50+4+10; %
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
        

    case '5Fork-Hairpin-70dT'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 4+70; %4dT spacer and 70dT tail;
        ConstructConst.dsDNA1 = 1710;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.ssDNA2 = 4+70+1719*2; %2*50+4+10; %
        ConstructConst.dsDNA2 = 1710;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 8]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
        
    case 'Hairpin-LambdaRH-dT-10-86bp'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 10; %0;
        ConstructConst.dsDNA1 = 3.265e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.265e3;
        ConstructConst.ssDNA2 = 2*86+4+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        %ConstructConst.fitrange = [4 12]; % Force in pN % set to [4 12] on 190308, had been [2 5]
        ConstructConst.fitrange = [5 10]; % Force in pN; set to [5 10] on 230816, though I'd been passing these values manually before. (AVT)
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
    
    case '5fork-10dT'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'Oligo Model';
        ConstructConst.dsDNA1 = 3431;
        ConstructConst.ssDNA1 = 4;
        ConstructConst.fitoffset = 1;
    %         fitrange = 200:260;
        ConstructConst.fitrange = [10 20]; % Force in pN
        ConstructConst.fakeoffset = 0;
        
        
    case 'T4'
        ConstructConst.hairpin = 0;
        %         dsDNA = 1.5e3;%first try
        %         dsDNA1 = 3011+29+120;%second try, ds part + overhang + stuck part
        ConstructConst.dsDNA1 = 3.4e3;
%         ConstructConst.dsDNA1 = 3.4e3 + 920;
%         ConstructConst.dsDNA1 = 1575;
%         ConstructConst.dsDNA1 = 4361-229+187;
        ConstructConst.ssDNA1 = 0;
        ConstructConst.fitoffset = 1;
        ConstructConst.fitrange = [4 7];
%         fitrange = [2 30]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 15];
        
    case 'Oligo'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'Oligo Model';
        ConstructConst.dsDNA1 = 3263;
        ConstructConst.ssDNA1 = 35;
        ConstructConst.fitoffset = 1;
%         fitrange = 200:260;
        ConstructConst.fitrange = [2 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        
    case 'fork'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'Oligo Model';
        ConstructConst.dsDNA1 = 3105;
        ConstructConst.ssDNA1 = 4;
        ConstructConst.fitoffset = 1;
%         fitrange = 200:260;
        ConstructConst.fitrange = [10 20]; % Force in pN
        ConstructConst.fakeoffset = 0;
        
     case 'RecA'
        ConstructConst.hairpin = 0;
        ConstructConst.dna1label = 'RecA Model';
        ConstructConst.dsDNA1 = 3263;
        ConstructConst.ssDNA1 = 140;
        ConstructConst.fitoffset = 1;
        ConstructConst.fitrange = [4 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        
    case 'Hairpin-RH-dT-19'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 19;
        ConstructConst.dsDNA1 = 3.05e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.05e3;
        ConstructConst.ssDNA2 = 19+89*2+4;
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitrange = [2 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [3 12];
        ConstructConst.stdF2 = [17.5 25];
        
    case 'Hairpin-RH-dT-38'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 38;
        ConstructConst.dsDNA1 = 3.05e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.05e3;
        ConstructConst.ssDNA2 = 2*89+4+38;%220
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitrange = [2 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case 'Hairpin-RH-dT-10'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 10; %0;
        ConstructConst.dsDNA1 = 3.055e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.055e3;
        ConstructConst.ssDNA2 = 2*89+4+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitrange = [5 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
        
    case 'Hairpin-2x-RH-dT-10'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 20; %0;
        ConstructConst.dsDNA1 = 3.055e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.055e3;
        ConstructConst.ssDNA2 = 2*89+4+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitrange = [5 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case 'Hairpin-LambdaRH-dT-10'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 10; %0;
        ConstructConst.dsDNA1 = 3.25e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 3.25e3;
        ConstructConst.ssDNA2 = 2*89+4+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 10]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case 'nick_test'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 1010; %0;
        ConstructConst.dsDNA1 = 2.25e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 2.25e3;
        ConstructConst.ssDNA2 = 2*86+4+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [2 5]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
    case 'GQ'
        ConstructConst.hairpin = 1;
        ConstructConst.dna1label = 'Hairpin Closed Model';
        ConstructConst.ssDNA1 = 34; %0;
        ConstructConst.dsDNA1 = 2.26e3;
        ConstructConst.dna2label = 'Hairpin Open Model';
        ConstructConst.dsDNA2 = 2.26e3;
        ConstructConst.ssDNA2 = 22+ConstructConst.ssDNA1; %2*50+4+10; %
        ConstructConst.fitoffset = 1;
%         fitrange = 100:150;
        ConstructConst.fitclosed = 1;
        ConstructConst.fitrange = [3 6]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 10];
        ConstructConst.stdF2 = [15 16];
        
     case 'MarcodsDNA'
        ConstructConst.hairpin = 0;
        %         dsDNA = 1.5e3;%first try
        ConstructConst.dsDNA1 = 1.565e3;
        ConstructConst.ssDNA1 = 0;
        ConstructConst.fitoffset = 1;
%         fitrange = 200:350;
        ConstructConst.fitrange = [2 30]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 15];
        
    case 'Saurabh_ssDNA'
        ConstructConst.hairpin = 0;
        %         dsDNA = 1.5e3;%first try
        ConstructConst.dsDNA1 = 0;
        ConstructConst.ssDNA1 = 1181;
        ConstructConst.fitoffset = 1;
        %         fitrange = 200:350;
        ConstructConst.fitrange = [2 30]; % Force in pN
        ConstructConst.fakeoffset = 0;
        ConstructConst.stdF = [2 15];
        
end