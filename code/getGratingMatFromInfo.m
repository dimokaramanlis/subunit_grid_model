function [stimmat] = getGratingMatFromInfo(stiminfo, spX, spY)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Nstimuli  = size(stiminfo, 1);
stimmat   = zeros(numel(spY) * numel(spX), Nstimuli, 'single');

[XX, YY] = meshgrid(spX, spY);
for istim = 1:Nstimuli
    spf     = stiminfo(istim,1);
    ortheta = stiminfo(istim,2);
    phangle = stiminfo(istim,3);
    stimmat(:,istim) = sin(2*pi*spf*(XX(:)*cos(ortheta) + YY(:)*sin(ortheta)) + phangle);
end

stimmat  = reshape(stimmat, numel(spY), numel(spX), Nstimuli);

end
