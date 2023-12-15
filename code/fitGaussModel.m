function mdlparams = fitGaussModel(mdlparams, stiminfo, stimorder, yy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
mdlfit = mdlparams;
rmax   = mdlparams.rmax/2;
%==========================================================================
% setup guess from current params
guess = [mdlfit.gaussparams mdlfit.outparams double(mdlfit.ktwts)];
%==========================================================================
% paramerer bounds 
lb = [guess(1) - rmax  guess(2) - rmax       1      1 -pi/4 -Inf   0 ...
    -Inf(1, numel(mdlparams.ktwts))];
ub = [guess(1) + rmax  guess(2) + rmax  rmax/2 rmax/2  pi/4  Inf Inf ...
    Inf(1, numel(mdlparams.ktwts))];
%==========================================================================
% actual optimization
options = optimoptions('fmincon','Algorithm','interior-point',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients', false,...
    'HessianFcn', [], 'FiniteDifferenceType', 'central', ...
     'MaxIterations',100);

foptim         = @(p) optimizeGaussianModel(p, stiminfo, stimorder, yy, mdlfit);

[fitprms, res] = fmincon(foptim, guess, [], [], [], [], lb, ub, [], options);
%==========================================================================
% return new structure
mdlparams.gaussparams = fitprms(1:5);
mdlparams.outparams   = fitprms(6:7);
mdlparams.ktwts       = fitprms(8:end);
%==========================================================================
end


function [f, g] = optimizeGaussianModel(gmparams, xx, xorder, yy, otherparams)

params              = otherparams;
params.gaussparams  = gmparams(1:5); 
params.outparams    = gmparams(6:7); 
params.ktwts        = gmparams(8:end); 

[lf, lg] = funfitGaussModel(params, xx, xorder);    

Np       = sum(yy);

f = -(log(lf)'*yy - sum(lf))/Np; %objective f

if nargout > 1
    g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
end


end
