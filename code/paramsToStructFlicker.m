function cellmdl = paramsToStructFlicker(mdlparams, gridcent, ktbas, dgrid, multfac)
%PARAMS_TO_STRUCT_FLICKER Converts model parameters to a structured format for flicker stimuli.
%
%   This function takes a vector of model parameters (mdlparams), 
%   information about the subunit grid (gridcent, dgrid), temporal 
%   basis functions (ktbas), and a scaling factor (multfac) to create a 
%   structured array (cellmdl) containing the model parameters in a more 
%   organized format. Subunits are sorted by distance to the grid center.
%
%   Inputs:
%       mdlparams - A vector of model parameters.
%       gridcent  - The center coordinates of the subunit grid [x, y].
%       ktbas     - Matrix of temporal basis functions.
%       dgrid     - The distance between subunits in pixels.
%       multfac   - A scaling factor relative to 7.5 um.
%
%   Outputs:
%       cellmdl   - A structured array containing the model parameters:
%                       cellmdl.subcnts   : Subunit centers [x, y] for each subunit
%                       cellmdl.gridcent  : Center of the grid [x, y]
%                       cellmdl.subwts    : Subunit weights
%                       cellmdl.ktwts     : Weights for temporal basis functions
%                       cellmdl.subsurrwt : Surround weights for temporal basis functions
%                       cellmdl.subsigma  : Subunit sigma (standard deviation)
%                       cellmdl.subsurrsc : Subunit surround scaling
%                       cellmdl.subbase   : Subunit nonlinearity midpoint
%                       cellmdl.outparams : Output nonlinearity parameters
%                       cellmdl.rmax      : Maximum response (80 * multfac)
%                       cellmdl.ktbasis   : Temporal basis functions (ktbas)
%                       cellmdl.kt        : Temporal kernel for subunit center (ktbas * cellmdl.ktwts)


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

