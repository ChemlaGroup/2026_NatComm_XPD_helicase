% 230816, AVT
% After return_FEC calculates extension and force, can optionally plot with
% plot_precalculated_FEC.
% Both programs derived from plot_offset_joint. This code plots fewer
% things.


function plot_precalculated_FEC(DNAext, FM, FA, FB)

global construct
global WLC_param

if nargin <3
   FA=[];
   FB=[];
end

ConstructConst = ConstructConstants(construct);

set(0,'defaultfigureposition',[100  200  1050  659]') % Put figures in middle of laptop screen.
max_force = max(FM)+1;
F = 0:0.05:max_force;
% [WLC_param] = WLC_parameters;
bp = 1*WLC_param.hds*XWLCContour(F,WLC_param.Pds,WLC_param.Sds,WLC_param.kT); % nm/bp
nt = 1*WLC_param.hss*XWLCContour(F,WLC_param.Pss,WLC_param.Sss,WLC_param.kT); % nm/nt


if ConstructConst.hairpin == 1
    % closed model
    WLC_ext = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp;
    % open model
    WLC_ext_op = ConstructConst.ssDNA2.*nt + ConstructConst.dsDNA2.*bp;
    figure; hold on;
    plot(WLC_ext, F, '--r', 'LineWidth', 2, 'DisplayName', 'WLC closed')
    plot(WLC_ext_op, F, '--k', 'LineWidth', 2, 'DisplayName', 'WLC open')
else
    WLC_ext = ConstructConst.ssDNA1.*nt + ConstructConst.dsDNA1.*bp;
    figure; hold on;
    plot(WLC_ext, F, '--r', 'LineWidth', 2, 'DisplayName', 'WLC model')
end

plot(DNAext, FB, 'Color', [0.6 0.6 0.75], 'LineWidth', 1, 'DisplayName', 'Trap B force')
plot(DNAext, FA, 'Color', [0.35 0.7 0.95], 'LineWidth', 1, 'DisplayName', 'Trap A force')
plot(DNAext, FM, 'b', 'LineWidth', 1, 'DisplayName', 'Average force')
xlabel('Extension (nm)');
ylabel('Force (pN)');
legend('location','northwest')
xlim([800 1200])
ylim([0 20])

end