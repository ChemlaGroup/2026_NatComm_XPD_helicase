% AVT, 230829
% Sub-routine for part of force_offset_subtraction: 
% offset sorting and smoothing.
% Saving as fields the sorted and final smoothed offset.

% 230830: instead of movmean, use smooth function, with option to use
% different methods.

% sample usage:
% offsetData=offset_smoothing(offsetData, trapData.sampperiod, 'movmean', 11);
% offsetData=offset_smoothing(offsetData, trapData.sampperiod);

function offsetData=offset_smoothing(offsetData, sampperiod, method, window)

if nargin <3
    method = 'moving'; %movmean
    window = [];
end

% Sort offset
    [offsetData.s_extOffset, sid] = sort(offsetData.TrapSep);
    offsetData.s_AOffsetX = offsetData.A_X(sid);
    offsetData. s_BOffsetX = offsetData.B_X(sid);
    offsetData.s_AOffsetY = offsetData.A_Y(sid);
    offsetData.s_BOffsetY = offsetData.B_Y(sid);

    

    % Determine window for smoothing
    if isempty(window)
        
        %Nfilt = round(1 / (trapData.sampperiod * 100)); % Divide by 10 to avoid oversmoothing
        % Above line actually does not filter, as written. This was apparently changed by MM from
        % earlier code, RemoveOffsetJoint, where filtering was taking place.
       
        Nfilt = round(1 / (sampperiod * 10)); % Divide by 10 to avoid oversmoothing
        %Nfilt = round(1 / (sampperiod * 30));
        
        % round to nearest odd number, for movmean to work properly
        Nfilt2 = 2*floor((Nfilt/2))+1;
        if ~isequal(mod(Nfilt2, 2), 1)
            warning('Offset filter window not odd; may cause artefacts!');
        end
    else
        Nfilt2=window;
    end
    
    %disp(Nfilt2)
    
    % Smooth offset
    if strcmp(method, 'median')
        ss_extOffset = movmedian(offsetData.s_extOffset,  Nfilt2);
        ss_AOffsetX = movmedian(offsetData.s_AOffsetX, Nfilt2);
        ss_BOffsetX = movmedian(offsetData.s_BOffsetX, Nfilt2);
        ss_AOffsetY = movmedian(offsetData.s_AOffsetY, Nfilt2);
        ss_BOffsetY = movmedian(offsetData.s_BOffsetY, Nfilt2);
    else
        ss_extOffset = smooth(offsetData.s_extOffset,  Nfilt2, method);
        ss_AOffsetX = smooth(offsetData.s_AOffsetX, Nfilt2,method);
        ss_BOffsetX = smooth(offsetData.s_BOffsetX, Nfilt2, method);
        ss_AOffsetY = smooth(offsetData.s_AOffsetY, Nfilt2, method);
        ss_BOffsetY = smooth(offsetData.s_BOffsetY, Nfilt2, method);
    end
    % My earlier version: movmean
%     ss_extOffset = movmean(offsetData.s_extOffset, Nfilt2);
%     ss_AOffsetX = movmean(offsetData.s_AOffsetX, Nfilt2);
%     ss_BOffsetX = movmean(offsetData.s_BOffsetX, Nfilt2);
%     ss_AOffsetY = movmean(offsetData.s_AOffsetY, Nfilt2);
%     ss_BOffsetY = movmean(offsetData.s_BOffsetY, Nfilt2);
    
    Nfilt3=round(Nfilt2/2);
    
    % Discard first halfwindow points (not as necessary for movmean as filter, but keep)
    sss_extOffset = ss_extOffset(Nfilt3:end);
    sss_AOffsetX = ss_AOffsetX(Nfilt3:end);
    sss_BOffsetX = ss_BOffsetX(Nfilt3:end);
    sss_AOffsetY = ss_AOffsetY(Nfilt3:end);
    sss_BOffsetY = ss_BOffsetY(Nfilt3:end);
    
    [offsetData.ssss_extOffset, index] = unique(sss_extOffset);
    offsetData.sss_AOffsetX = sss_AOffsetX(index);
    offsetData.sss_BOffsetX = sss_BOffsetX(index);   
    offsetData.sss_AOffsetY = sss_AOffsetY(index);
    offsetData.sss_BOffsetY = sss_BOffsetY(index);

end