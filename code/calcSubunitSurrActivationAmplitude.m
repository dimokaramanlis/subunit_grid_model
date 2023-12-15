function [Fc, Fs, Jc, Js] = calcSubunitSurrActivationAmplitude(params, spf)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
spff = spf(:); %standardize input
sigma = params.subsigma;
kscale = params.subsurrsc;
%==========================================================================
Fc = exp(-(pi * sqrt(2) * sigma * spff).^2);
Fs = Fc.^(kscale^2);
%==========================================================================
if nargout > 2
    Jc = - 4 * pi^2  * sigma   * Fc .* (spff.^2);
    J1 = - 4 * pi^2 * kscale^2 * sigma   * Fs .* (spff.^2);
    J2 = - 4 * pi^2 * kscale   * sigma^2 * Fs .* (spff.^2);
    Js = [J1 J2];
end

%==========================================================================
end

