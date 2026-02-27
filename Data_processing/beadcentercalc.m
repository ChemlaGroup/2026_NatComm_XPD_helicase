function [results] = beadcentercalc(im, pltimgyn, saveimgyn, filename, figdir)

if nargin < 2
pltimgyn = 0 ; 
end

if nargin < 3
saveimgyn = 0 ; 
end

if nargin < 4
filename = 'bead image'; 
end

titlename = filename; 

% figure();imshow(im); 
% figure();imshow(imsub); 
% figure();imshow(imsub1); 
% figure();imshow(imsub2); 

im_sub = im;
rowsum = sum(im_sub,2); colsum = sum(im_sub,1); 
rx = 1:length(rowsum); cx = 1:length(colsum); 

% figure()
% subplot(5,5,[1:4])
% plot(cx,colsum,'k'); hold on; 
% xlim([1 length(colsum)]);set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
% subplot(5,5,[10:5:25])
% plot(rx,rowsum); hold on; 
% xlim([1 length(rowsum)]); view(90,90);set(gca,'YTickLabel',[]);;set(gca,'XTickLabel',[]);
% subplot(5,5,[6:9 11:14 16:19 21:24])
% imshow(imsub2); 
          
[val,idx] = max( im_sub(:));
[y,x] = ind2sub( [size(im_sub)], idx ); 

approxres.cmu = x; approxres.rmu = y; 
approxres.cs = 100; approxres.rs = 100; pixmin = 25;
[xp,yp] = meshgrid(approxres.cmu-min(pixmin, approxres.cs) :approxres.cmu+min(pixmin, approxres.cs) ,...
            approxres.rmu-min(pixmin, approxres.rs) :approxres.rmu+min(pixmin, approxres.rs) );
im_crop = imcrop(im_sub, ...
            [approxres.cmu-min(pixmin, approxres.cs), approxres.rmu-min(pixmin, approxres.rs), 2*min(pixmin, approxres.cs), 2*min(pixmin, approxres.rs)]);
[results] = autoGaussianSurf(xp,yp,im_crop); 

% figure()
% imshow(im_crop)

imc = im_crop;
rowsum = sum(imc,2); colsum = sum(imc,1); 
rmin = approxres.rmu-min(pixmin, approxres.rs); rmax = approxres.rmu+min(pixmin, approxres.rs); 
rx = rmin:rmax; 
cmin = approxres.cmu-min(pixmin, approxres.cs); cmax = approxres.cmu+min(pixmin, approxres.cs);
cx = cmin:cmax; 

results.pcorr = corr(imc(:), results.G(:)); 

mu = results.x0 ; sig = results.sigmax ; vo = 0 ; 
% amp = results.a;
amp = max(colsum);
cfx = cmin:0.1:cmax;
cf = amp*gaussmf(cfx,[sig,mu]); 

mu = results.y0 ; sig = results.sigmay ; vo = 0 ; 
% amp = results.b; 
amp = max(rowsum);
rfx = rmin:0.1:rmax;
rf = amp*gaussmf(rfx,[sig,mu]);


    if pltimgyn == 1
    % close all
    figure();
    
    ax1 = subplot(5,5,[1:4]);
    plot(cx,colsum,'k'); hold on; plot(cfx,cf,'r');
    xlim([cmin cmax])
    set(gca,'YTickLabel',[]);set(gca,'xaxisLocation','top')
    set(gca,'XTickLabel',[]);
    
    ax2 = subplot(5,5,[10:5:25]);
    plot(rx,rowsum); hold on; plot(rfx,rf,'r');
    xlim([rmin rmax])
    view(90,90);
    set(gca,'YTickLabel',[]);set(gca,'xaxisLocation','top')
    set(gca,'XTickLabel',[]);
    
    ax3 = subplot(5,5,[6:9 11:14 16:19 21:24]);
    imshow(imc); hold on; 
    title(titlename ,'Interpreter','none')
    
    
        if saveimgyn == 1 
        figname = [ titlename  '_imgfit'];
        saveas(gcf,[figdir figname '.fig']); 
        saveas(gcf,[figdir figname '.png']);
        else
            return
        end
    
    
    figure()
    surf(imc); hold on; 
    surf(results.G);
    view(2);
    axis off
    title(titlename ,'Interpreter','none')
    
    
        if saveimgyn == 1 
        figname = [ titlename '_imgcompare'];
        saveas(gcf,[figdir figname '.fig']); 
        saveas(gcf,[figdir figname '.png']);
        else
            return
        end
   
    % pause
    % close all

    else 
        return

    end


res=[results.x0 results.y0];

end
