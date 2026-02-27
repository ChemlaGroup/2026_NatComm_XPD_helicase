
% Find maximum of each burst for each trace in 'data' and save in
% burst2peak field. 

function data = burst_peak(data)

for i=1:length(data)    % process each trace separately
    Nbursts = size(data{i}.burst2, 2);
    data{i}.burst2peak = [];
    
    for j = 1:Nbursts      % for each burst in a trace, find maximum point and record the corresponding time
        ind = find(data{i}.time > data{i}.burst2(1, j) & data{i}.time < data{i}.burst2(2, j));
        [bpmax, indmax] = max(data{i}.bp(ind)); 
        %[bpmax, indmax] = max(movingmeanAT(data{i}.bp(ind), 11)); 
        data{i}.burst2peak(j) = data{i}.time(ind(indmax));
    end

end


end

