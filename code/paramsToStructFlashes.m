function cellmdl = paramsToStructFlashes(mdlparams, gridcent, dgrid)
%UNTITLED Summary of this function goes here
%   dgrid is distance between subunits in pixels
 
Nsubs     = numel(mdlparams) - 9; 
mainpts  = generateHexSubunitGrid(Nsubs);
 
cellmdl.subcnts   = gridcent + mainpts*dgrid;
cellmdl.gridcent  = gridcent;
cellmdl.subwts    = mdlparams(1:Nsubs)';
cellmdl.subsigma  = mdlparams(Nsubs + 1);
cellmdl.subsurrsc = mdlparams(Nsubs + 2);
cellmdl.subsurrwt = mdlparams(Nsubs + 3);
cellmdl.subparams = mdlparams(Nsubs + (4:5));
cellmdl.outparams = mdlparams(Nsubs + (6:8));
cellmdl.outbase   = mdlparams(Nsubs + 9);
cellmdl.cwindow   = 60;
        
end

