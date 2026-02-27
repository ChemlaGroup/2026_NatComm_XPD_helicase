% 241003, AVT

% Making this the new version of getBetaFactors.
% Instead of reading in beta taken on same day as measurements, return beta
% from the day that the calibration grid was imaged, for data after
% steerable mirror adjustment.
% Leaving old processing as before for now (using daily beta).


function beta = getBetaFactors(startpath, Date)

%no_mixing = 0;

numDate = str2double(Date);

if numDate > 240909
    imDate = '240920';
    pathname = 'enter_correct_filepath';
    fixedxyres = load([pathname 'fixedxyimdir_0/fixedxyres_' imDate '.mat']);
    %fixedxyres = load([pathname 'fixedxyimdir_1/fixedxyres_' imDate '.mat']);
    beta = fixedxyres.beta(1:2, 1:2);
    
    % Leave old method untouched for now
else 
    imDate = Date;
    
    try
        fixedxyres = load([startpath 'fixedxyres_' Date '.mat']);
    catch
        fixedxyres=[];
        disp('No beta file to load.')
    end
            
    if ~isempty(fixedxyres)
        beta = fixedxyres.beta(1:2, 1:2);
        
        %Look for bead images to calculate beta directory/subdirectories
    else
        files = dir([startpath '**/*GridMove*.bmp']);
        
        if ~isempty(files)
            [output, ~, ~] = findfixedxy(Date);
            beta = output.beta;
        else
            % if no fixedxyres can be determined, use beta with no mixing
            beta =   eye(2); %dummy identity matrix
            disp('Cannot access or calculate beta factors. Assuming no mixing.')
        end
        
    end
    
end

disp(['Using beta from ' num2str(imDate) ', if available.'])
disp(beta)

end