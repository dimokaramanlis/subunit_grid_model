function [fpred, insignal] = predictSubGridSubSurrModel(mdlparams, stiminfo, stimorder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

kt     =  mdlparams.ktbasis *  mdlparams.ktwts;
ktsurr =  mdlparams.ktbasis *  mdlparams.subsurrwt;

Nt = numel(kt);
Nstim = size(stimorder, 2);
%Nsubs = numel(mdlparams.subwts);


iuse = mdlparams.subwts>0;
Nsubs = nnz(iuse);

ArgSubs    = cos(calcSubunitActivationPhases(mdlparams, stiminfo));
[Rfc, Rfs] = calcSubunitSurrActivationAmplitude(mdlparams, stiminfo(:, 1));
subcentacts  = Rfc .* ArgSubs(:,iuse);
subsurracts  = Rfs .* ArgSubs(:,iuse);

tempcentacts  = reshape(subcentacts(stimorder, :), [Nt, Nstim * Nsubs]);
tempsurracts  = reshape(subsurracts(stimorder, :), [Nt, Nstim * Nsubs]);

%tempfin  = kt' * tempcentacts;
tempfin  = kt' * tempcentacts + ktsurr' * tempsurracts;
tempfin  = reshape(tempfin, [size(stimorder,2), Nsubs]);

RsubsNlin    = reshape(logisticfun(tempfin + mdlparams.subbase),[Nstim, Nsubs]);
% insignal     = RsubsNlin * mdlparams.subwts + mdlparams.outparams(1);
% fpred        = mdlparams.outparams(2) * logisticfun(insignal);% get spiking response
insignal     = RsubsNlin * mdlparams.subwts(iuse);
% outparams = [mdlparams.outparams;1];
% fpred = softplusfun3m(outparams, insignal);

fpred = nakarushton(mdlparams.outparams, insignal);

fpred     = gather(fpred);

end

