function  Rsqf = sufTrialRsq( trialCounts, prs)
% sufTrialRsq computes a measure of trial-to-trial consistency (R-squared) 
% across different stimulus combinations.
%
% Inputs:
%   trialCounts - A matrix of spike counts, where each row corresponds to a 
%                 cell and each column corresponds to a stimulus presentation.
%                 Size: [Ncells x Nstim].
%   prs         - A vector indicating the stimulus combination for each column 
%                 in `trialCounts`. Each unique value corresponds to a different 
%                 stimulus combination.
%
% Outputs:
%   Rsqf - A vector of averaged R-squared values (one for each cell) that quantifies 
%          the consistency of responses across even and odd trials for each 
%          stimulus combination.
%
% Function Details:
%   - The function partitions the spike count data into odd and even trials for
%     each stimulus combination.
%   - It computes the R-squared values comparing the means of even and odd 
%     trials, assessing the trial-to-trial consistency.
%   - Two R-squared metrics are calculated:
%       Rsq: Compares odd trial means to even trial means.
%       Rsq2: Compares even trial means to odd trial means.
%   - The final R-squared (`Rsqf`) is the average of Rsq and Rsq2 for each cell.
%
% Usage:
%   This function is useful for evaluating the consistency of neural responses 
%   across repeated presentations of a stimulus in a given experimental setup.

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

