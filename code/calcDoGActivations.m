function [Fc, Fs, Jc, Js] = calcDoGActivations(params, xvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
spfreq = xvals(:, 1); 
ori    = xvals(:, 2);
%==========================================================================
[Rfc,   Rjc] = calcGaussianActivationAmplitude(params.gaussparams, [spfreq ori]);
[Rfs,   Rjs] = calcGaussianSurrActivationAmplitude(params, [spfreq ori]);
[argf, argj] = calcGaussianActivationPhase(params.gaussparams, xvals);
Fc           = Rfc .* cos(argf);
Fs           = Rfs .* cos(argf);
%==========================================================================
if nargout > 2
    
    Jc = zeros(size(xvals,1), 6);
    Js = zeros(size(xvals,1), 6);
    
    Jc(:, 1:2) =  -Rfc .* sin(argf) .* argj;
    Jc(:, 3:5) = cos(argf) .* Rjc;
    
    Js(:, 1:2) =  -Rfs .* sin(argf) .* argj;
    Js(:, 3:6) = cos(argf) .* Rjs;

end
%==========================================================================
end

