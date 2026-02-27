% 250609, AT

% Downsample fleezers raw data to as near old trap sampling frequency (100 Hz) as
% possible.
% Intended to be run on data immediately after loading (via ReadJointFile)
% before subtraction of offset, calculation of force, extension, etc.

function data = fl_filter_decimate(data)

ds_fields = {'A_X', 'A_Y', 'A_Sum', 'B_X', 'B_Y', 'B_Sum', 'A_FB_X', 'A_FB_Y', 'A_FB_Sum',...
    'B_FB_X', 'B_FB_Y', 'B_FB_Sum', 'C_X', 'C_Y', 'C_Y_Sum', 'C_X_Sum', 'time', 'trappos1', 'trappos2'};

dfactor = ceil(0.01/data.sampperiod);


for j=1:numel(ds_fields)
    
    if isfield(data, ds_fields{j})
        data.(ds_fields{j}) = FilterAndDecimate(data.(ds_fields{j}), dfactor);
    else
        disp(['Field ' ds_fields{j} ' not found.']);
    end    
    
end

data.sampperiod=data.sampperiod*dfactor;
data.downsampfactor = dfactor;

end