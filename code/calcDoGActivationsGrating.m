function [Fout, J] = calcDoGActivationsGrating(params, xvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
spfreq = xvals(:, 1); 
ori    = xvals(:, 2);
surrwt    = params.surrwt;
%==========================================================================
[Rfc,   Rjc] = calcGaussianActivationAmplitude(params.gaussparams, [spfreq ori]);
[Rfs,   Rjs] = calcGaussianSurrActivationAmplitude(params, [spfreq ori]);
[argf, argj] = calcGaussianActivationPhase(params.gaussparams, xvals);
Fout         = (Rfc - surrwt * Rfs) .* cos(argf);
%==========================================================================
if nargout > 1
    
    J12  = (- Rfc + surrwt * Rfs) .* sin(argf) .* argj;
    J345 = cos(argf) .* (Rjc - surrwt * Rjs(:, 1:3));
    J6   = cos(argf) .* (- surrwt * Rjs(:, 4));
    J7   = -Rfs .* cos(argf);
    J = [J12 J345 J6 J7];
end

%==========================================================================
end

