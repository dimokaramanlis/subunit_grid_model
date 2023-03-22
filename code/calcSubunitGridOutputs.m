function [F,  J] = calcSubunitGridOutputs(params, xstimuli)
%CALCSUBUNITGRIDGRATINGRESPONSE Summary of this function goes here
%   Detailed explanation goes here
%   xstimuli : N x 4 matrix where N is number of stimuli with 4 parameters
%               width, orientation, phase, contrast bias
% TO GET THE SUBUNIT WEIGHTS, THE WEIGHT GRID WITH ITS RESPECTIVE LOCATIONS
% MUST BE PROVIDED

%SUBCENTERS COME AS (X, Y) pairs from (0, 0) center
%TURN THEM TO REAL CENTERS first 

%make params a structure, and unpack it

%==========================================================================
spfreq   = xstimuli(:, 1); 
[Rf, Rj] = calcSubunitActivationAmplitudes(params, spfreq);
argf     = calcSubunitActivationPhases(params, xstimuli);
F        = Rf .* cos(argf);
%==========================================================================
if nargout > 1
    J = reshape(Rj, [size(xstimuli,1), 1, size(Rj, 2)]) .* cos(argf);
end

%==========================================================================
end

