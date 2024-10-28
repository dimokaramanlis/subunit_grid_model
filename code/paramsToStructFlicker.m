function cellmdl = paramsToStructFlicker(mdlparams, gridcent, ktbas, dgrid, multfac)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 % calculate activations for every regularization
 

Nwt       = size(ktbas,2);
Nsubs     = numel(mdlparams) - 6 - 2*Nwt; 
mainpts   = generateHexSubunitGrid(Nsubs);
 
cellmdl.subcnts = gridcent + mainpts*dgrid;
cellmdl.gridcent  = gridcent;

cellmdl.subwts    = mdlparams(1:Nsubs)';
cellmdl.ktwts     = mdlparams(Nsubs + (1:Nwt))';
cellmdl.subsurrwt = mdlparams(Nsubs + Nwt + (1:Nwt))';
cellmdl.subsigma  = mdlparams(Nsubs + 2 * Nwt + 1);
cellmdl.subsurrsc = mdlparams(Nsubs + 2 * Nwt + 2);
cellmdl.subbase   = mdlparams(Nsubs + 2 * Nwt + 3);
cellmdl.outparams = mdlparams(Nsubs + 2 * Nwt + (4:6))';
cellmdl.rmax      = 80 * multfac;
cellmdl.ktbasis   = ktbas;
cellmdl.kt        = ktbas * cellmdl.ktwts;


end

