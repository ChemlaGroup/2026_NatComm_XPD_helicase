% Originally by Monika Makurath, based in part on previous lab codes. --AVT

% modified by AVT (August 2023):
% -- 2D offset subtraction: subtract offset both in X and in Y regardless of dominant pulling direction.
% Note: no longer have separate DNAext_X and DNAext_Y fields; single
% DNAext field.

% -- re-introduced smoothing (factor had been set to 1 in
%earlier code). Replaced filter function with movmean, using odd
%window.

function trapData = force_offset_subtraction(offsetData, trapData, calData, param, method, window)

if nargin <5
    method = 'moving'; %movmean
    window = [];
end


% QPD OFFSET SUBTRACTION
trapData.A_X = trapData.A_X-calData.AXoffset; % both are in volts
trapData.B_X = trapData.B_X-calData.BXoffset; % both are in volts
trapData.A_Y = trapData.A_Y-calData.AYoffset; % both are in volts
trapData.B_Y = trapData.B_Y-calData.BYoffset; % both are in volts
offsetData.A_X = offsetData.A_X-calData.AXoffset; % both are in volts
offsetData.B_X = offsetData.B_X-calData.BXoffset; % both are in volts
offsetData.A_Y = offsetData.A_Y-calData.AYoffset; % both are in volts
offsetData.B_Y = offsetData.B_Y-calData.BYoffset; % both are in volts


if param.offset_file == 1
    % Sort and smooth offset
    offsetData=offset_smoothing(offsetData, trapData.sampperiod, method, window);
    
    % Interpolate offset values corresponding to data values
    ssi_AOffsetX = interp1(offsetData.ssss_extOffset, offsetData.sss_AOffsetX, trapData.TrapSep);
    ssi_BOffsetX = interp1(offsetData.ssss_extOffset, offsetData.sss_BOffsetX, trapData.TrapSep);
    ssi_AOffsetY = interp1(offsetData.ssss_extOffset, offsetData.sss_AOffsetY, trapData.TrapSep);
    ssi_BOffsetY = interp1(offsetData.ssss_extOffset, offsetData.sss_BOffsetY, trapData.TrapSep);
    
    % Replace the raw values with force offset subtracted values   
    trapData.A_X = trapData.A_X - ssi_AOffsetX;
    trapData.B_X = trapData.B_X - ssi_BOffsetX;
    trapData.A_Y = trapData.A_Y - ssi_AOffsetY;
    trapData.B_Y = trapData.B_Y - ssi_BOffsetY;
    
end

%% CONVERT QPD DATA FROM VOLTS TO NANOMETERS
trapData.A_Xnm = (trapData.A_X)*calData.alphaAX;  % convert QPD data from volts to nm
trapData.B_Xnm = (-trapData.B_X)*calData.alphaBX; % convert QPD data from volts to nm
trapData.A_Ynm = (-trapData.A_Y)*calData.alphaAY; % convert QPD data from volts to nm
trapData.B_Ynm = (trapData.B_Y)*calData.alphaBY;  % convert QPD data from volts to nm

%% DNA EXTENSION IN NANOMETERS

BeadSep_X = trapData.TrapSep_X - trapData.A_Xnm - trapData.B_Xnm;
BeadSep_Y = trapData.TrapSep_Y - trapData.A_Ynm - trapData.B_Ynm;
trapData.DNAext = sqrt(BeadSep_X.^2 + BeadSep_Y.^2) - calData.beadA/2 - calData.beadB/2;

%% FORCE (pN)
trapData.force_AX = trapData.A_Xnm * calData.kappaAX; % multiply QPD data by trap stiffness
trapData.force_BX = trapData.B_Xnm * calData.kappaBX; % multiply QPD data by trap stiffness
trapData.force_AY = trapData.A_Ynm * calData.kappaAY; % multiply QPD data by trap stiffness
trapData.force_BY = trapData.B_Ynm * calData.kappaBY; % multiply QPD data by trap stiffness

end


