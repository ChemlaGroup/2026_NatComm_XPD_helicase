% Calculates bp unwound for hairpin construct. Calls XWLCCountor.m

% AVT: I got this version from Barbara.
% Note: only accurate if fakeoffset has been taken into account.

function nbp = bpunwound(F,X,Favg,params)


global WLC_param

if nargin <4
    %params is an array of the following:
    lds = 3050; %length of ds handles
    lss = 10; %length of ss loading site
    lhp = 89; %length of hairpin
    lloop = 4; %length of T loop in hairpin
else
    lds = params(1);
    lss = params(2);
    lhp = params(3);
    lloop = params(4);
end
if nargin <3
    Favg = nanmean(F);
end
    
%lhp %To check that construct info is being processed correctly 
%ds DNA WLC params
Pds = WLC_param.Pds;
Sds = WLC_param.Sds;
hds = WLC_param.hds;

%ss DNA WLC params
Pss =WLC_param.Pss;
Sss = WLC_param.Sss;
hss = WLC_param.hss;


hairpinwidth = 2;
ssclosed = lss;
ssopen = lss+2*lhp+lloop;


xclosed = lds*hds*XWLCContour(F,Pds,Sds) + ssclosed*hss*XWLCContour(F,Pss,Sss) + hairpinwidth;
xopen = lds*hds*XWLCContour(F,Pds,Sds) + ssopen*hss*XWLCContour(Favg,Pss,Sss); %length of DNA at mean force with hairpin fully open

idx1 = find(X < xopen);
idx2 = find(X >= xopen);
nbp = zeros(1,length(X));

if (~isempty(idx2)) % if hairpin opens all the way
    
    nbp(idx1) = (X(idx1) - xclosed(idx1))./(2*hss*XWLCContour(F(idx1),Pss,Sss)); %teth_alph.fit(3)
    
    nbp(idx2) = (X(idx2) - xclosed(idx2))./(2*hss*XWLCContour(F(idx2),Pss,Sss));
    
else
    
    nbp = (X - xclosed)./(2*hss*XWLCContour(F,Pss,Sss));
    
end

end