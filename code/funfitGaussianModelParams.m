function [F, J] = funfitGaussianModelParams(params, xvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
outparams = params.outparams;
%==========================================================================
% get inner prediction
[fin, jin] = calcGaussianActivationsGrating(params.gaussparams, xvals);

% get spiking response
[F, jout]   = rlogistic3(params.outparams, fin);
%==========================================================================
if nargout > 1
    prefact = outparams(2)*F.*(1-F/outparams(3));

    J1to5  = prefact .* jin;
    J      = [J1to5 jout];
end
%==========================================================================
end

