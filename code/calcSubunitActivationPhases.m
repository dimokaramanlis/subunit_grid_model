function [F] = calcSubunitActivationPhases(params, spforiphase)
%The function can deal with mulitple subunits
%   Detailed explanation goes here

%Returns derivatives for grid spacing, xshift and yshift
%==========================================================================
spf = spforiphase(:, 1); ori = spforiphase(:, 2); sph = spforiphase(:, 3);
%==========================================================================
xo   = params.subcnts(:,1);
yo   = params.subcnts(:,2);

fact1   = sqrt(xo.^2 + yo.^2)';

inangle = ori - atan(yo./xo)';
fact2   = cos(inangle);

F = (2*pi*spf) .* (fact1 .* fact2) + (sph - pi/2);
%==========================================================================
end

