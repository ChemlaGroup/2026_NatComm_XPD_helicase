% Function from Steve Yeo

% Modified:
% 230905: replacing slashes with global fbslash (AVT)
% 230908: read image files from subdirectory if they have been moved (AVT)

% 230910: Only offer option to move files if they are still in data
% directory (AVT). Output saved to directory with data.

function [output, res, results] = findfixedxy(date)

global datadirectories
global fbslash

date = num2str(date);
strdate = date; 
datadir = [datadirectories date fbslash];

% If files have already been moved to a subdirectory, check latest, then first:
sub_imdir = dir([datadir 'fixedxyimdir*']);
if ~isempty(sub_imdir)
     if ~isempty(dir([datadir sub_imdir(end).name fbslash '*GridMove*']))
        imdir = [datadir sub_imdir(end).name fbslash];
     elseif ~isempty(dir([datadir sub_imdir(1).name fbslash '*GridMove*']))
        imdir = [datadir sub_imdir(1).name fbslash];
     else
        imdir = datadir;
    end
else
    imdir = datadir;
end
disp(imdir);
extra_title = ' fixed x/y based on grid positions';
imagesavename = [date '_fixedxy'];


%% get background 
imamp = 1; % 1 or 255
clear x ; x =[]; 
files = dir([imdir '*background*.bmp']); 
for i = 1:length(files)  

    im = imread([imdir files(i).name]);
    im = imamp*im2double(im); 
    x(:,:,i) = im; 
end
bgp = mean(x,3); 
disp('background done'); 
%% get moving images 
xarr = [];
yarr = []; 
files = dir([imdir '*GridMove*.bmp']); 
for i = 1:length(files) 
        strarr = split(files(i).name, '_');
        xarr = [xarr, strarr(3)]; 
        yarr = [yarr, strarr(4)]; 
end

[xi,yi] = meshgrid(1:640,1:480);
xarr = unique(xarr) ;
yarr = unique(yarr) ;

%% get fixed trap bead position 
disp('fixed trap calc.. ')
clear x ; x = []; 
files = dir([imdir '*fixed*.bmp']); 
for i = 1:length(files) 
    disp(i); 
    im = imread([imdir files(i).name]);
    im = imamp*im2double(im); 
    x(:,:,i) = im; 
%     figure()
%     imshow(im); 
%     pause
%     close all 
end
medim = mean(x,3);

imsub = imsubtract(medim, bgp); 
imsub = imsub./(max(imsub(:))); % normalize
imsub1 = imsub.*(imsub>0); % only above 0 data 
thr = mean(imsub1(:))+2*std(imsub1(:)) ;  % threshold 
imsub2 = imsub.*(imsub>thr); % only above threshold 
% get approximate center 
%try
[fixed_results] = beadcentercalc(imsub2,1,0); 
%end
%% get movable trap bead positions 
disp('movable trap calc.. ')
res=[]; 
clear totimg; 
totimg = []; 
n = 1; 
for i = 1:length(xarr)
    for j = 1:length(yarr)
        files = dir([imdir '*_' char(xarr(i)) '_' char(yarr(j)) '_*.bmp']);
%        disp(['doing...... X: ' char(xarr(i)) ' , Y: ' char(yarr(j))])
        for k = 1:length(files)
%             disp(files(k).name); 
            im = imread([imdir files(k).name]);
            im = imamp*im2double(im); 
            x(:,:,k) = im; 
        end
        medim  = mean(x,3); 
        imsub = imsubtract(medim, bgp); % subtract background
        imsub = imsub./(max(imsub(:))); % normalize
        imsub1 = imsub.*(imsub>0); % only above 0 data 
        thr = mean(imsub1(:))+2*std(imsub1(:)) ;  % threshold 
        imsub2 = imsub.*(imsub>thr); % only above threshold 
        totimg(:,:,n) = imsub2;
        n=n+1; 
        % get approximate center 
        try
        [results] = beadcentercalc(imsub2,0); 
        res=[res; str2num(xarr{i}), str2num(yarr{j}), results.x0, results.y0, results.pcorr, results.ssep, results.r2];
        end
    end
end
%% show sum image 

% figure 
% medmvim = mean(totimg,3);
% summvim = sum(totimg,3);
% imshow(summvim)

%% Plot in volt 
xpos = res(:,3); 
idxx = find(xpos >= (mean(xpos)-2*std(xpos)) & xpos <= (mean(xpos)+2*std(xpos)) );
ypos = res(:,4); 
idxy = find(ypos >= (mean(ypos)-2*std(ypos)) & ypos <= (mean(ypos)+2*std(ypos)) );
iid = intersect(idxx, idxy);

% figure 
% scatter(res(iid,1)./100,res(iid,2)./100, 100, abs(res(iid,6)), 'filled') ; hold on; 
% xlabel('X(V)')
% ylabel('Y(V)')

%%
xpos = res(:,3); 
idxx = find(xpos >= (mean(xpos)-2*std(xpos)) & xpos <= (mean(xpos)+2*std(xpos)) );
ypos = res(:,4); 
idxy = find(ypos >= (mean(ypos)-2*std(ypos)) & ypos <= (mean(ypos)+2*std(ypos)) );
iid = intersect(idxx, idxy);
% iid = 1:length(xpos); 



figure
scatter(res(iid,3),res(iid,4), 100, abs(res(iid,6)), 'filled') ; hold on; 
c = colorbar; 
c.Label.String = 'sum squared error';
grid on; 
scatter(fixed_results.x0, fixed_results.y0, 500, 'rx'); 
xlabel('X pixels')
ylabel('Y pixels')
xlim([min(res(iid,3))-1 max(res(iid,3))+15]);
ylim([min(res(iid,4))-1 max(res(iid,4))+1]);
fixed_x = griddata(res(iid,3), res(iid,4), res(iid,1), fixed_results.x0, fixed_results.y0,'cubic');
fixed_y = griddata(res(iid,3), res(iid,4), res(iid,2), fixed_results.x0, fixed_results.y0,'cubic');
%clc; %230910: turning off clc, with fitting output suppressed (AVT) 
disp(fixed_x/100); 
disp(fixed_y/100); 
disp(sum(res(:,6)));

g = zeros(2, 1);
% g(1) = plot(NaN,NaN,'r--', 'displayname', '1 RPA model');
g(1) = plot(NaN,NaN,'ko','displayname', 'movable trap pos. ');hold on; 
g(2) = plot(NaN,NaN,'rx','displayname', 'fixed trap pos. ');hold off; 
legend(g); 
set(gca, 'YDir','reverse')
axis equal
title({[strdate ' ' extra_title] ['fixed x(V): ' num2str(fixed_x/100,'%01.6f') ' fixed y(V): ' num2str(fixed_y/100,'%01.6f')]})
%toc    

saveas(gcf, [datadirectories date fbslash imagesavename '.fig']);

output.movable = res; 
output.fixed = [fixed_results.x0, fixed_results.y0]; 
output.fixedV = [fixed_x/100, fixed_y/100]; 



%% get beta that converts V to pixel. Y = bX formulation. 

sortedpos  = sortrows(res, [ 1 2 ] ); 

VM = sortedpos(:,[1 2]); 
VM1 = [VM./100, ones(size(VM,1),1)]; % include offset. Since the values gotten above is for V/100 values get correct V values.  
PM = sortedpos(:,[3 4]);

beta = PM'/VM1'; % Get linear transformation 

output.beta = beta; 

%% save results 

assignin('base',"fixedxyres",output)
var = [ 'fixedxyres1' ];
eval([var  '= output' ]);   
save([datadir fbslash 'fixedxyres_' strdate '.mat'],'-struct', var)

%%
% Optionally move images to subfolder. If files already moved, skip.
if strcmp(imdir, datadir)
    prompts = {'Move Fixed images to new folder? (y or n)'};
    defaults = {'y'};
    conditions = inputdlg(prompts,['Move fixed x y images: '],1,defaults);
    disp(conditions);
    if conditions{1} == 'y'
        i=0
        % savedirbase = [imdir 'fixedxyimdir_' num2str(i) '/'];
        % imdir = [datadirectory date ];
        saveimgdir = ['fixedxyimdir_' num2str(i) fbslash];
        while exist([imdir saveimgdir], 'dir');
            i=i+1;
            saveimgdir = ['fixedxyimdir_' num2str(i) fbslash];
        end
        mkdir([imdir saveimgdir]);
        movefile([imdir '*' imagesavename '.fig'],[imdir saveimgdir],'f');
        movefile([imdir '*.bmp'],[imdir saveimgdir],'f');
    else
        disp('did not move fixed xy images and figures. ')
    end
    
end







end