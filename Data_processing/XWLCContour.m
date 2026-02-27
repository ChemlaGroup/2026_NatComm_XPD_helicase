function x = XWLCContour(F, P, S, kT, option)
% x = XWLCContour(F, P, S, kT, option)
%
% Calculates extension over contour x/L given force F (pN), 
% persistence length P (nm) and stretch modulus S (pN).
%
% If option = 1 returns extension over contour x/L;
% if option = 2 returns stiffness times contour KL.

% Default Values
if nargin < 4
    kT = 4.14;
end
if nargin < 3
    S = 1200;
end
if nargin < 2
    P = 50;
end

% Solves XWLC by making the simplification x/l = 1+F/S-1/u,
% where u^3+a*u+b = 0, and solves the cubic equation

% Simplification Variables
a = -4*(F*P/kT-0.75);
b = -4;
q = a/3; 
r = -b/2;
D = q.^3+r^2;

sgn = 2*(a <= 0)-1; %accounts for change of sign and uses correct soln (YRC 2/08)
% if a <= 0 sgn = 1, a > 0 sgn = -1;
u = (r+sqrt(D)).^(1/3)+sgn.*(sgn.*(r-sqrt(D))).^(1/3);
x = 1+F./S-1./u;

if nargin == 5 & option >= 1 %return dKL/dx
    w = 0.5./(1-x+F/S).^3+1;
    k = kT/P*w./(1+kT/P/S*w);
    dwdx = 1.5./(1-x+F/S).^4.*(1-k/S);
    dkdx = kT/P*dwdx./(1+kT/P/S*w).^2;
    if option == 1
        x = k;
    end;
    if option == 2
        x = dkdx;
    end;    
end; 