function [ Rsqf] = sufTrialRsq( trialCounts, prs)
%IMAGESELECTIVITYINDEX Calculates the selectivity to natural images
%   First defined in Quian Quiroga et al. (2007), also used by the Allen
%   Institue Brain Observatory.
%   Inputs:
%           allTrialCounts: n x Nt x p matrix, n is number of cells, Nt is number of trials
%   and p number of images
%   Outputs:
%           sIndex: n x 1 vector

[Ncells, Nstim] = size(trialCounts);
Ncombis         = max(prs);
%==========================================================================

evenMean = NaN(Ncells, Ncombis); 
oddMean  = NaN(Ncells, Ncombis); 

for icombi = 1:Ncombis
    comcounts = trialCounts(:, prs == icombi);
    oddMean(:, icombi)  = mean(comcounts(:, 1:2:end), 2);
    evenMean(:, icombi) = mean(comcounts(:, 2:2:end), 2);
end

Rsq  = 1 - sum((oddMean-evenMean).^2, 2) ./sum((oddMean  - mean(oddMean, 2)).^2,2);
Rsq2 = 1 - sum((oddMean-evenMean).^2, 2) ./sum((evenMean - mean(evenMean,2)).^2,2);
Rsqf = mean([Rsq Rsq2],2);
%==========================================================================
end

