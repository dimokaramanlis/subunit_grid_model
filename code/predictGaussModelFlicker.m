function fpred = predictGaussModelFlicker(mdlparams, stiminfo, stimorder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


rfactiv     = calcGaussianActivationsGrating(mdlparams.gaussparams, stiminfo);
allspactivs = rfactiv(stimorder);
allactivs   =  mdlparams.ktwts * (mdlparams.ktbasis' * allspactivs);


alpha = mdlparams.outparams(2);
gamma = mdlparams.outparams(1);
fpred = alpha * (1 + exp(-allactivs' - gamma)).^-1;


end

