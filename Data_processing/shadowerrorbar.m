% from Yann, 240205
% plots shaded bands in place of error bars

% AVT: made plots optional

function [px,py] = shadowerrorbar(x,y,y_errU,y_errL)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Remove any NaN's
idx_y = find(isnan(y));
idx_y_errU = find(isnan(y_errU));
idx_y_errL = find(isnan(y_errL));
idx = setxor(1:length(y),[idx_y idx_y_errU idx_y_errL]);

% Define perimeter of object to fill
py = [y(idx)+y_errU(idx) fliplr(y(idx)-y_errL(idx))];
px = [x(idx) fliplr(x(idx))];

makeplot=0;

if makeplot
    % Plot mean and fill area defined by error bars
    plot(x(idx),y(idx),'k'); hold on;
    fill(px, py, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'None');
end

end

