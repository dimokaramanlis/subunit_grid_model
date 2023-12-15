function mdlparams = fitDoGModel(mdlparams, stiminfo, stimorder, yy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
mdlfit = mdlparams;
rmax   = mdlparams.rmax/2;
mdlfit.surrktwts = -mdlfit.ktwts/20; %small surround strength
mdlfit.surrsc    = 2;
Nwts = numel(mdlparams.ktwts);
%==========================================================================
% setup guess from current params
guess = [mdlfit.gaussparams mdlfit.surrsc mdlfit.outparams ...
    double(mdlfit.ktwts) double(mdlfit.surrktwts)];
%==========================================================================
% paramerer bounds 
lb = [guess(1) - rmax  guess(2) - rmax       1      1 -pi/4  1 -Inf   0 ...
    -Inf(1, 2 * Nwts)];
ub = [guess(1) + rmax  guess(2) + rmax  rmax/2 rmax/2  pi/4 10 Inf Inf ...
     Inf(1, 2 * Nwts)];
%==========================================================================
% actual optimization
options = optimoptions('fmincon','Algorithm','interior-point',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients', false,...
    'HessianFcn', [], 'FiniteDifferenceType', 'central', ...
    'MaxIterations', 100);

foptim   = @(p) optimizeGaussianModel(p, stiminfo, stimorder, yy, mdlfit);
[fitprms, res] = fmincon(foptim, guess, [], [], [], [], lb, ub, [], options);
%==========================================================================
% return new structure
mdlparams.gaussparams = fitprms(1:5);
mdlparams.surrsc      = fitprms(6);
mdlparams.outparams   = fitprms(7:8);
mdlparams.ktwts       = fitprms(8 + (1:Nwts));
mdlparams.surrktwts   = fitprms(8 + Nwts + (1:Nwts));
%==========================================================================
end


function [f, g] = optimizeGaussianModel(gmparams, xx, xorder, yy, otherparams)


params              = otherparams;
Nwts                = numel(params.ktwts);

params.gaussparams  = gmparams(1:5); 
params.surrsc       = gmparams(6); 
params.outparams    = gmparams(7:8); 
params.ktwts        = gmparams(8 + (1:Nwts)); 
params.surrktwts    = gmparams(8 + Nwts + (1:Nwts)); 


[lf, lg] = funfitDoGModel(params, xx, xorder);    

Np       = sum(yy);

f = -(log(lf)'*yy - sum(lf))/Np; %objective f

if nargout > 1
    g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
end


end

function [c, ceq, gradc, gradceq] = filternlndog(p)


Nwts = (numel(p) - 8)/2;
ktc  = p(8 + (1:Nwts));
kts  = p(8 + Nwts + (1:Nwts));
c(1) = ktc*kts';
ceq = [];

if nargout > 2
    gradc = zeros(size(p));
    gradc(8 + (1:Nwts)) = kts;
    gradc(8 + Nwts + (1:Nwts)) = ktc;
    gradceq = [];
end


end

