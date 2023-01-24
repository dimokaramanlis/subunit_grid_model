function [F, J] = calcGaussianActivationAmplitude(gparams, spftheta)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
spff   = spftheta(:, 1); %standardize input
thori  = spftheta(:, 2); %standardize input

sigma1  = gparams(3);
sigma2  = gparams(4);
theta0  = gparams(5);
%==========================================================================
F = exp(- 2 * (pi * spff).^2 .* ...
    (sigma2^2 * sin(thori  + theta0).^2 + sigma1^2 * cos(thori + theta0).^2));
%==========================================================================
if nargout > 1
    % FIX!!!
    J1 = -2 * (pi * spff).^2 .* F .* (2 * sigma1 * cos(thori + theta0).^2);
    J2 = -2 * (pi * spff).^2 .* F .* (2 * sigma2 * sin(thori + theta0).^2);
    
    J3 = -2 * (pi * spff).^2 .* F .* sin(2*(thori + theta0)) * (sigma2^2 - sigma1^2);
    
    J  = [J1 J2 J3];
end
%==========================================================================
end

