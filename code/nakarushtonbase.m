function [ fout, J] = nakarushtonbase( params, x )
%LOGISTICFUN evaluates modified logistic with params at x
%   Detailed explanation goes here
%==========================================================================
%define inputs
x       = x(:);
Rmax    = params(1);
n       = params(2);
kappa   = params(3); 
b       = params(4);
%==========================================================================
%calculate function value
xp     = x.^n;
f      = Rmax * xp./(xp + kappa^n);
fout   = f + b;
%==========================================================================
if nargout>1
    J1 = f/Rmax;
    J2 = f .* (log(x) -  (log(x) .* xp + log(kappa) * kappa^n)./(xp + kappa^n));
    J3 = - f * n * kappa^(n-1) ./ (xp + kappa^n);
    J = [J1 J2 J3 ones(size(J1))];
end
%==========================================================================
end

