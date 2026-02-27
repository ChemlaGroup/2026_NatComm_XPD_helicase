% from Steve Yeo; copied Sept. 2024
% updated

function ConvFactor = getConvFactor(StrDate)
% disp('getconvfactor')

numDate = str2double(StrDate); 


if numDate < 240410
ConvFactor.flconv = 123;
elseif numDate > 240410
 ConvFactor.flconv= 123; % new conv factor 
end 



if numDate < 180629
ConvFactor.otconvX =  891.5;
ConvFactor.otconvY = 505.7;
elseif numDate < 230810
ConvFactor.otconvX =  883.59;
ConvFactor.otconvY = 508.93;
elseif numDate < 240909      %numDate  > 230810
 ConvFactor.otconvX = 896.6691;
ConvFactor.otconvY = 534.2339;
elseif numDate > 240909 % Use temporarily while OT's steerable mirror out of sorts
    ConvFactor.otconvX = 862.8293; %888.6375;
    ConvFactor.otconvY=0; % steerable mirror not functioning in Y
end 











end