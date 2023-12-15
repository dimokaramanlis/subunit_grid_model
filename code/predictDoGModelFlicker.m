function fpred = predictDoGModelFlicker(mdlparams, stiminfo, stimorder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



[fc, fs]  = calcDoGActivations(mdlparams, stiminfo);

allfc     = fc(stimorder);
allfs     = fs(stimorder);
cbaseproj = (mdlparams.ktbasis' * allfc);
sbaseproj = (mdlparams.ktbasis' * allfs);

allactiv    = mdlparams.ktwts * cbaseproj +  mdlparams.surrktwts * sbaseproj;
allactiv    = allactiv';



alpha = mdlparams.outparams(2);
gamma = mdlparams.outparams(1);
fpred = alpha *(1 + exp(-(allactiv + gamma))).^-1;

end

