% from BKS

% part of basicanalysis4alice

% Pick out start and end locations of bursts

function data=burst_selection_BKS(data)

for k = 1:length(data)  %length(data)
    num = 2; %number of subplots
     %if ~isfield(data{k},'burst2')
    data{k}.burst2 = [];
    
    figure; ax = subplot(num,2,[1,2]); %plot(data{k}.time,data{k}.bp,'k');
    hold on; 
    %plot(downsample(data{k}.time,10),downsample(data{k}.bp,10),'k'); ylim([-5, 90])
plot(downsample(data{k}.time,2),downsample(data{k}.bp,2),'k'); ylim([-5, 90]) %changing 201104 for low-force data

    title([num2str(k) '  ' data{k}.date '\_' num2str(data{k}.refnums(4)) '/n Select Burst Start and End Times one burst at a time']); ylabel('Base Pairs Open');
    set(gca,'YTick',0:5:100,'ygrid','on');
    
%     ax = [ax subplot(num,2,num*2-1:num*2)]; plot(downsample(data{k}.time,10),downsample(data{k}.force,10),'b');
%     ylim([0, 15]); ylabel('Force (pN)'); set(gca,'ygrid','on');
%     
%     ax = [ax subplot(num,2,num*2-1:num*2)]; plot(downsample(data{k}.apd.time,20),downsample(data{k}.apd.apd2,20),'Color',[0,0.5,0]);
%     ylabel('Cy3 Intensity'); set(gca,'ygrid','on');
%     
    ax = [ax subplot(num,2,num*2-1:num*2)]; plot(data{k}.time,data{k}.rawdata.data.A_Y,'r'); hold all; plot(data{k}.time,data{k}.rawdata.data.B_Y);
    
    linkaxes(ax,'x');
    xlabel('Time (s)');
    
    s = 0;
    while s == 0
        pause;
        [x,y] = ginput(2);
        
        if length(x) == 2 % if only time point is selected, try again
            if x(1)>x(2) % if second time point is earlier than first, clear all bursts
                data{k}.burst2 = [];
            else % if 2 chronological time points are selected, add to list of bursts
                data{k}.burst2 = [data{k}.burst2, x];
                plot(data{k}.burst2(:,end),[1;1],'bv');
            end
        elseif isempty(x) % if nothing selected move on to next trace
            s = 1;
        end
    end
    
    close(gcf);
    data{k}.proc2 = [];
    
    for i = 1:size(data{k}.burst2,2)
        ind = find(data{k}.time >= data{k}.burst2(1,i) & data{k}.time <= data{k}.burst2(2,i));
        
        if ~isempty(nonzeros((data{k}.bp(ind)>25)))
            data{k}.proc2 = [data{k}.proc2, 1];
        else
            data{k}.proc2 = [data{k}.proc2, 0];
        end
     end
    %end
end
end
