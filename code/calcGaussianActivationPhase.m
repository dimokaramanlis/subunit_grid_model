function [F, J] = calcGaussianActivationPhase(gparams, spforiphase)
%The function can deal with mulitple subunits
%   Detailed explanation goes here

%Returns derivatives for grid spacing, xshift and yshift
%==========================================================================
spf = spforiphase(:, 1); 
ori = spforiphase(:, 2); 
sph = spforiphase(:, 3);
%==========================================================================
xo      = gparams(1);
yo      = gparams(2);
%==========================================================================
fact1   = sqrt(xo^2 + yo^2);
inangle = ori - atan(yo/xo);
fact2   = cos(inangle);

F       = 2*pi*spf .* fact1 .* fact2 + sph - pi/2;
%==========================================================================
if nargout > 1
    
    sinfact = sin(inangle);
    
    df1xo = xo/fact1;
    df2xo = -sinfact * yo/fact1.^2;
    J1    = 2*pi*spf .* (df1xo * fact2 + fact1 * df2xo);
    
    df1yo = yo/fact1;
    df2yo = sinfact * xo/fact1.^2;
    J2    = 2*pi*spf .* (df1yo * fact2 + fact1 * df2yo);
    
    J = [J1 J2];
end
%==========================================================================
end

