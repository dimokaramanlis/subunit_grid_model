function mdlparams = gfFitDoGMonotonic(mdlparams, xvals, yy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
mdlfit = mdlparams;
rmax   = mdlparams.rmax/2;
%==========================================================================
% setup guess from current params
guess = [mdlfit.gaussparams mdlfit.surrsc mdlfit.surrwt mdlfit.outparams];
%==========================================================================
% paramerer bounds 
lb = [guess(1) - rmax  guess(2) - rmax       1      1 -pi/4 ...
    1   0 -Inf -Inf  0];
ub = [guess(1) + rmax  guess(2) + rmax  rmax/2 rmax/2  pi/4  ...
    6   1 Inf  Inf Inf];
%==========================================================================
% actual optimization
options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients', false,...
    'HessianFcn', [], 'FiniteDifferenceType', 'central');

foptim         = @(p) optimizeDoGModel(p, xvals, yy, mdlfit);

[fitprms, res] = fmincon(foptim, guess, [], [], [], [], lb, ub, [], options);
%==========================================================================
% return new structure

mdlparams.gaussparams = fitprms(1:5);
mdlparams.surrsc      = fitprms(6);
mdlparams.surrwt      = fitprms(7);
mdlparams.outparams   = fitprms(8:end);

%==========================================================================
end


function [f, g, H] = optimizeDoGModel(gmparams, xx, yy, otherparams)

params             = otherparams;
params.gaussparams = gmparams(1:5); 
params.surrsc      = gmparams(6);
params.surrwt      = gmparams(7);
params.outparams   = gmparams(8:end);

[lf, lg] = gfmodels.funfitDoGModelParams(params, xx);    
%[lf] = gfmodels.funfitDoGModelParams(params, xx);    

Np       = sum(yy);

f = -(log(lf)'*yy - sum(lf))/Np; %objective f

if nargout > 1
    g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
end

if nargout > 2
    ht1 = ((yy./(lf.^2)).*lg)' * lg;
    ht2 = ((yy./lf)' * lg + lg' * (yy./lf))/sum(wtparams);
    ht3 = -(sum(lg, 1) + sum(lg, 1)')/sum(wtparams);
    H = (ht1 + ht2 + ht3) / Np;
end


end
