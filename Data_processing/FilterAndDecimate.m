% Shared by Yann, 250529

function x = FilterAndDecimate(xd, num)

x=filter(ones(1, num)/num, 1, xd);
idx = num:num:length(xd);

x=x(idx);