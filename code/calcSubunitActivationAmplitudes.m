function [F, J] = calcSubunitActivationAmplitudes(params, spf)
%The function can deal with mulitple subunits
%   Detailed explanation goes here

%Returns derivatives for grid spacing, xshift and yshift
%==========================================================================
spff = spf(:); %standardize input
sigma  = params.subsigma;
surrsc = params.subsurrsc;
surrwt = params.subsurrwt;
%==========================================================================

Fc = exp(-(pi * sqrt(2) * sigma)^2 * spff.^2);

Fs = Fc.^(surrsc^2);
%F  = Fc + surrwt * Fs;
F  = Fc - surrwt * Fs;

%==========================================================================
if nargout > 1
    
    J1 = -4 * pi^2 * sigma * (spff.^2) .* (Fc - surrwt * surrsc^2 * Fs);
    J2 = (surrwt * 4 * (pi^2) * (sigma^2) * surrsc) * Fs .* (spff.^2) ;
    J3 = - Fs;
    
    J = [J1 J2 J3];
    
end
%==========================================================================
end

