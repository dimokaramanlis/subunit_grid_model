function [Fout, J] = calcGaussianActivationsGrating(params, xvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
spfreq = xvals(:, 1); 
ori    = xvals(:, 2);
%==========================================================================
[Rf,    Rj]  = gfmodels.calcGaussianActivationAmplitude(params, [spfreq ori]);
[argf, argj] = gfmodels.calcGaussianActivationPhase(params, xvals);
Fout         = Rf .* cos(argf);
%==========================================================================
if nargout > 1
    
    J12  = (-Rf.* sin(argf)) .* argj;
    J345 = cos(argf) .* Rj;
    
    J = [J12 J345];
end

%==========================================================================
end

