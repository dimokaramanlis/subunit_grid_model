function generators = getGaussGeneratorsFixMovie(experiment, ...
            imEnsemble, runningfixations, modelparams, xmar, ymar)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nx  = experiment.projector.screen(2);
Ny  = experiment.projector.screen(1);
Nyx = Ny * Nx;
[Nyim, Nxim, Nimages]  = size(imEnsemble);
%imEnsemble = flip(imEnsemble, 1);
%runningfixations(2, :, :) = Nyim - runningfixations(2, :, :);
%--------------------------------------------------------------------------
Nblocks = size(runningfixations, 3);
Nframes = size(runningfixations, 2);
Nt      = size(modelparams.ktbasis, 1);
%--------------------------------------------------------------------------
cparams    = modelparams.gaussparams;
cparams(1) = cparams(1) + xmar + 0.5;
cparams(2) = Ny - (cparams(2) + ymar + 0.5);

% cparams(1) = cparams(1) + xmar ;
% cparams(2) = Ny - (cparams(2) + ymar);

cparams(5) = pi - cparams(5);

cx   = round(cparams(1));
cy   = round(cparams(2));
rmax = round(modelparams.rmax/2);
rx   = cx + (-rmax:rmax);
ry   = cy + (-rmax:rmax);

blockstimulus = zeros([numel(ry), numel(rx), Nframes], 'single');
generators    = zeros(Nframes - Nt + 1, Nblocks, 'single');

[xx, yy]   = meshgrid(rx, ry);

spcomp     = gauss2dfun([cparams 1], {xx(:), yy(:)});
spcomp     = spcomp/norm(spcomp);
kt         = modelparams.ktbasis * modelparams.ktwts;
kt         = kt/norm(kt);


for iblock = 1:Nblocks 
    
    fixrun = runningfixations(:, :, iblock);
    %------------------------------------------------------------------
    % generate stimulus
    for ifix = 1:size(fixrun, 2)
        [xmin, xmax, ymin, ymax, xOr, yOr] = getRanges(...
            fixrun(2, ifix), fixrun(3, ifix), Nx, Ny, Nxim, Nyim, rx, ry);
        blockstimulus(ymin:ymax, xmin:xmax, ifix) = imEnsemble(yOr, xOr, fixrun(1, ifix));
    end
    %------------------------------------------------------------------
    % reduce with spatial receptive field
    spsignals = spcomp' * reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 
    %------------------------------------------------------------------
    % convolve with temporal filter
    generators(:, iblock) = conv(spsignals, flip(kt), 'valid');
    %------------------------------------------------------------------
end
generators = reshape(generators, [(Nframes - Nt + 1) * Nblocks, 1]);

end

function [xmin, xmax, ymin, ymax, xOr, yOr] = getRanges(trX, trY, Nxs, Nys, Nx, Ny, rx, ry)
     
ymin = (Nys / 2) - trY + 1; 
if ymin <= 1, ymin = 1; end

ymax = Ny + (Nys / 2) - trY; 
if (ymax > Nys),  ymax = Nys; end

xmin = (Nxs / 2) - trX + 1; 
if xmin <= 1, xmin = 1; end

xmax = Nx + (Nxs / 2) - trX; 
if (xmax > Nxs),  xmax = Nxs; end

% Rx    = xmin:xmax;
% rxuse = ismembc(rx, Rx);
% 
% Ry    = ymin:ymax;
% ryuse = ismembc(ry, Ry);

rxuse = (rx >= xmin & rx <=xmax);
ryuse = (ry >= ymin & ry <=ymax);

xOr = trX - (Nxs / 2) + rx(rxuse);
yOr = trY - (Nys / 2) + ry(ryuse);

xmin = find(rxuse, 1, 'first');
xmax = find(rxuse, 1, 'last');
ymin = find(ryuse, 1, 'first');
ymax = find(ryuse, 1, 'last');

end

