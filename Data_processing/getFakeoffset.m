function fakeoffset = getFakeoffset(DNAext, FM, fitrange, fake_offset_0, construct)

global WLC_param
ConstructConst = ConstructConstants(construct);
% 190117 Steve 
% Get fakeoffset for FEC by : 
% 1) Get force array of data within fitrange 
% 2) Get ext from data corresponding to 1) 
% 3) Get ext from dmodel corresponding to 1)
% 4) Minimize the difference. 


fkoffsetidx = find(FM > fitrange(1) & FM < fitrange(2)); % get FM index within the fitrange. 
fkoffset_Farr = FM(fkoffsetidx); % get FM force array that's within the fitrange 
DNAext_offset = DNAext(fkoffsetidx) ; % get DNAext corresponding to the fkoffsetidx


% get WLC extension using fkoffset_Farr ( force array of the data. ) 
bp = 1*WLC_param.hds*XWLCContour(fkoffset_Farr, WLC_param.Pds, WLC_param.Sds, WLC_param.kT); % nm/bp
nt = 1*WLC_param.hss*XWLCContour(fkoffset_Farr, WLC_param.Pss, WLC_param.Sss, WLC_param.kT); % nm/nt
% closed model
WLC_ext_fit = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp; 
% open model
%     WLC_ext_op_fit = ConstructConst.ssDNA2.*nt + ConstructConst.dsDNA2.*bp;

% find o where the difference is minized. 
fakeoffset = fminsearch(@(x) sum((WLC_ext_fit - (DNAext_offset - x)).^2),fake_offset_0);
