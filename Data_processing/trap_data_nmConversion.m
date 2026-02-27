% Originally by Monika; updated by AVT

% 230802: adding TrapSep_X and TrapSep_Y as variables
% 230818: saving above as fields
% 230910: implementing beta and gamma conversion for OT.
% -- beta: 2*2 matrix of V-to-pixel conversion factors, based on bead grid images
% -- gamma: nm/pixel conversion factor from image of 10 um grid

%240916: 
% --adding Steve's implementation of beta (from his trap_data code)

%240918
% -- updating OTconv values with Steve's approach; pulling from file
% depending on date. 


function trapData = trap_data_nmConversion(trapData, param, beta)

if nargin <3
    beta = eye(2); % dummy identity matrix
end

ConvFactor = getConvFactor(param.Date);

%flconv = 123; 
%otconvX = 883.59; % updated on 06-29-2018
%otconvY = 508.93; % updated on 06-29-2018
%ot_gamma = 144.8; % nm/pixel, average for x and y (based on 2018 grid image analysis by SY)
flconv = ConvFactor.flconv; 
otconvX = ConvFactor.otconvX;
otconvY = ConvFactor.otconvY;


if param.instrument == 0
    %finish when gamma value for fleezers available
    trapData.TrapSep = (trapData.trappos2 - trapData.trappos1)*flconv;
    trapData.TrapSep_X=trapData.TrapSep;
    trapData.TrapSep_Y=trapData.TrapSep*0;
elseif param.instrument == 1
    
    %From Steve's code:
        
        pixperV_X = norm(beta*[1; 0]); pixperV_Y = norm(beta*[0; 1]); % get 1V step in x and y to get pixel per volts 
        nmperpix_X = otconvX/pixperV_X; nmperpix_Y = otconvY/pixperV_Y; % get nm/px. oTconv = nm/px( from cal grid ) * px/V ( from bead grid ).
        %disp(['nmperpix ' nmperpix_X ' in X, ' nmperpix_Y 'in Y'])
        TrapSep_Vx = trapData.trappos1 - trapData.t2x ; TrapSep_Vx = reshape(TrapSep_Vx, 1, []);
        TrapSep_Vy = trapData.trappos2 - trapData.t2y ; TrapSep_Vy = reshape(TrapSep_Vy, 1, []);
        TrapSep_V = [TrapSep_Vx; TrapSep_Vy]; % if sm tilted, sm V to pixel frame will include both x and y components 
        TrapSep_P = beta*TrapSep_V; % convert V to Pix using beta 
        
        %temporary work-around for OT data in 1D:
        if isnan(nmperpix_Y)
            nmperpix_Y=145.8333;
            nmperpix_X = 140.6250; % first peak value
             %nmperpix_X = 148.4375 % mean peak value --gives worse result
            
            disp(['Using nmperpix factors from ' param.Date '.'])
            %disp(['nmperpix_X: ' num2str(nmperpix_X)  '; otconvX: ' num2str(otconvX)])
        end
        
        trapData.TrapSep_X = TrapSep_P(1,:)*nmperpix_X; % Trapsep V to nm 
        trapData.TrapSep_Y = TrapSep_P(2,:)*nmperpix_Y; % Trapsep V to nm
        % Parameters were traceData.TrapSep_Xnm and traceData.TrapSep_Ynm
        % in Steve's code. Taking off 'nm' for compatibility with the rest of my code. 
        trapData.TrapSep = sqrt(trapData.TrapSep_X.^2+trapData.TrapSep_Y.^2); 
        % in Steve's code, what I call TrapSep above was TrapSepXY, and the
        % old TrapSep left for backward compatibility.
        %traceData.TrapSep = sqrt(((traceData.trappos1 - traceData.t2x)*otconvX).^2 + ((traceData.trappos2 - traceData.t2y)*otconvY).^2);

    
    %disp([gamma_X, gamma_Y]);
    %trapData.TrapSep_X = Trap_Sep_pixels(1, :)*gamma_X;
    %trapData.TrapSep_Y = Trap_Sep_pixels(2, :)*gamma_Y;
    %trapData.TrapSep = sqrt(trapData.TrapSep_X.^2 + trapData.TrapSep_Y.^2);
    
end