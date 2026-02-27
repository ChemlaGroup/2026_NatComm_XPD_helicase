function dG = SeqFreeEnergy2(seq,loop,model,salt)

% Free energies [kcal/mol] in order of (5'-XY-3'):
%        AA/TT AC/TG AG/TC AT/TA CA/GT CC/GG CG/GC GA/CT GC/CG TA/AT
dG_UO = [-1.27 -1.71 -1.53 -1.12 -1.72 -2.08 -2.50 -1.57 -2.53 -0.84];
dG_FR = [-1.23 -1.49 -1.36 -1.17 -1.66 -1.93 -2.37 -1.47 -2.36 -0.84];

% Salt correction factors [kcal/mol]
m_UO = 0.114*ones(1,10);
m_FR = [0.145 0.099 0.070 0.117 0.091 0.063 0.132 0.155 0.079 0.091];

% Remove spaces
idx = findstr(seq,' ');
seq = seq(setxor(1:length(seq),idx));

% Default model and salt concentration  
if nargin <= 3
    salt = 0.020 + 0.1 + 8.0*sqrt(0.003); % [M] Empirical formula: [Na+]eq = [Mon+] + 3.3*sqrt([Div++])
    if nargin == 2
        model = 'FR';
    end    
end

% Convert to free energies [kT] at specific salt concentration
mol = 602.2; % 1 mol = 602.2e21
kcal = 4184; % 1 kcal = 4184 J
kT = 1.381e-2*298; % 1 kT = 1.38e-23*298 J 
switch model 
    case 'FR'
        dG_NN = (dG_FR - m_FR*log(salt))*1.712; %kcal/mol/kT;
    case 'UO'
        dG_NN = (dG_UO - m_UO*log(salt))*1.712; %kcal/mol/kT;
end

% Nearest Neighbor free energies
for i = 1:length(seq)-1
    twobase = seq(i:i+1);
    % AA/TT AC/TG AG/TC AT/TA CA/GT CC/GG CG/GC GA/CT GC/CG TA/AT
    switch twobase
        case {'AA','TT'}
            dG(i) = dG_NN(1);
        case {'AC','GT'}
            dG(i) = dG_NN(2);
        case {'AG','CT'}
            dG(i) = dG_NN(3);
        case 'AT'
            dG(i) = dG_NN(4);
        case {'CA','TG'}
            dG(i) = dG_NN(5);
        case {'CC','GG'}
            dG(i) = dG_NN(6);
        case 'CG'
            dG(i) = dG_NN(7);
        case {'GA','TC'}
            dG(i) = dG_NN(8);
        case 'GC'
            dG(i) = dG_NN(9);
        case 'TA'
            dG(i) = dG_NN(10);
    end
end    
       
if length(loop) > 0
    dG(length(seq)) = 2.43*1.712; %kcal/mol/kT;
    dG((length(seq)+1):(length(seq)+length(loop))) = 0;
end    
