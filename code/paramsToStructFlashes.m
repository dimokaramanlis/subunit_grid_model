function cellmdl = paramsToStructFlashes(mdlparams, gridcent, dgrid)
%PARAMS_TO_STRUCT_FLASHES Converts model parameters to a structured format.
%
%   This function takes a vector of model parameters (mdlparams) and 
%   additional information about the subunit grid (gridcent, dgrid) to 
%   create a structured array (cellmdl) containing the model parameters 
%   in a more organized and accessible format. Subunits are sorted by 
%   distance to the grid center.
%
%   Inputs:
%       mdlparams - A vector of model parameters. The first Nsubs elements
%                   are the subunit weights.
%       gridcent  - The center coordinates of the subunit grid [x, y].
%       dgrid     - The distance between subunits in pixels.
%
%   Outputs:
%       cellmdl   - A structured array containing the model parameters:
%                     cellmdl.subcnts   : Subunit centers [x, y] for each subunit
%                     cellmdl.gridcent  : Center of the grid [x, y]
%                     cellmdl.subwts    : Subunit weights
%                     cellmdl.subsigma  : Subunit sigma (standard deviation)
%                     cellmdl.subsurrsc : Subunit surround scaling
%                     cellmdl.subsurrwt : Subunit surround weight
%                     cellmdl.subparams : Subunit nonlinearity parameters 
%                     cellmdl.outparams : Output nonlinearity parameters (a, n, k)
%                     cellmdl.outbase   : Output nonlinearity baseline(b)
%                     cellmdl.cwindow   : Rectangual window (set to 60, not used) 
%
 
Nsubs    = numel(mdlparams) - 9; 
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

