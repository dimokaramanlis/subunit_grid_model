function [stimmat, stimparams] = getGratingFlashStimulus(stimdesc, spX, spY, bwflag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Nwidths   = numel(stimdesc.stripewidths);
Nstimuli  = sum(stimdesc.Nphases .* stimdesc.Norientations);

stimparams = zeros(Nstimuli, 3);
stimmat    = zeros(numel(spY) * numel(spX), Nstimuli, 'uint8');

[XX, YY] = meshgrid(spX, spY);

istim = 1;

for iw = 1:Nwidths
    
    spf   = 0.5/stimdesc.stripewidths(iw);
    nph   = stimdesc.Nphases(iw);
    nori  = stimdesc.Norientations(iw);
    
    for iori = 1:nori
        
        ortheta = pi*(iori-1)/nori;
        
        
        for iphase = 1:nph
            
            phangle = 2*pi*(iphase-1)/nph;
            vals = sin(2*pi*spf*(XX(:)*cos(ortheta) + YY(:)*sin(ortheta)) + phangle);
            if bwflag; vals = sign(vals); end
            stimmat(:, istim) = uint8(255 * (vals + 1)/2);
            stimparams(istim, :) = [spf ortheta phangle];
            istim = istim + 1;
            
        end 
    end
end

stimmat  = reshape(stimmat, numel(spY), numel(spX), Nstimuli);

end
