function generators = generatorsFixationMovie(modelparams, ...
            imEnsemble, runningfixations, xmar, ymar)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
% let's figure out the model type

%--------------------------------------------------------------------------

Nx  = modelparams.screenx;
Ny  = modelparams.screeny;
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
cparams(5) = pi - cparams(5);

cx   = round(cparams(1));
cy   = round(cparams(2));
rmax = round(modelparams.rmax/2);
rx   = cx + (-rmax:rmax);
ry   = cy + (-rmax:rmax);

blockstimulus = zeros([numel(ry), numel(rx), Nframes], 'single');
generators    = zeros(Nframes - Nt + 1, Nblocks, 'single');

[xx, yy]   = meshgrid(rx, ry);

% case Gauss:
spcomp     = gauss2dfun([cparams 1], {xx(:), yy(:)});
spcomp     = spcomp/norm(spcomp);
ktwts      =  modelparams.ktwts(:);
kt         = modelparams.ktbasis *ktwts;
kt         = kt/norm(kt);

% case DoG

spcompcent = gauss2dfun([cparams 1], {xx(:), yy(:)});
spnorm     = norm(spcompcent);
spcompcent = spcompcent/spnorm;

surrparms = cparams;
surrparms(3:4) = surrparms(3:4) * modelparams.surrsc;
spcompsurr = gauss2dfun([surrparms 1], {xx(:), yy(:)});
spcompsurr = spcompsurr/spnorm/ modelparams.surrsc^2;

kt         =  modelparams.ktbasis * modelparams.ktwts;
ktsurr     = modelparams.ktbasis *  modelparams.surrktwts;
cnorm      = norm(kt);
kt         = kt/cnorm;
ktsurr     = ktsurr/cnorm;

% subunit grid


for ii = 1:Nmodels
    modelparams(ii).subcnts(:, 1) = modelparams(ii).subcnts(:, 1) + xmar + 0.5;
    modelparams(ii).subcnts(:, 2) = Ny - (modelparams(ii).subcnts(:, 2) + ymar + 0.5);
    
    
    
    modelparams(ii).kt     = gpuArray(modelparams(ii).ktwts'     *...
        modelparams(ii).ktbasis');
    modelparams(ii).ktsurr = gpuArray(modelparams(ii).subsurrwt' *...
        modelparams(ii).ktbasis');
    
    pxarr = round(200*modelparams(ii).subsigma);
    ikeep = modelparams(ii).subwts > 0;

    xcents = modelparams(ii).subcnts(ikeep,1);
    ycents = modelparams(ii).subcnts(ikeep,2);
    subfilters = (xx(:) - xcents').^2 + (yy(:) - ycents').^2;
    subfilters(sqrt(subfilters) > pxarr) = NaN; % filter useless values

    subfilterscent = exp(-subfilters/(2* modelparams(ii).subsigma^2));
    subfilterscent(isnan(subfilters)) = 0;
    subfilterscent   = subfilterscent./(2 * pi * modelparams(ii).subsigma^2);

    subfilterssurr = exp(-subfilters/...
        (2* (modelparams(ii).subsurrsc * modelparams(ii).subsigma)^2));
    subfilterssurr(isnan(subfilterssurr)) = 0;
    subfilterssurr   = subfilterssurr./...
        (2 * pi *(modelparams(ii).subsurrsc * modelparams(ii).subsigma)^2);

    modelparams(ii).subfilterscent = gpuArray(subfilterscent);
    modelparams(ii).subfilterssurr = gpuArray(subfilterssurr);
    
    sumall = sum(modelparams(ii).subwts);
    subwts = modelparams(ii).subwts(ikeep)*sumall/sum(modelparams(ii).subwts(ikeep));
    modelparams(ii).subwts = gpuArray(subwts);
end



for iblock = 1:Nblocks 
    
    fixrun = runningfixations(:, :, iblock);
    %------------------------------------------------------------------
    % generate stimulus
    for ifix = 1:size(fixrun, 2)
        [xmin, xmax, ymin, ymax, xOr, yOr] = getRanges(...
            fixrun(2, ifix), fixrun(3, ifix), Nx, Ny, Nxim, Nyim, rx, ry);
        blockstimulus(ymin:ymax, xmin:xmax, ifix) = imEnsemble(yOr, xOr, fixrun(1, ifix));
    end
    %======================================================================
    % gaussian
    % reduce with spatial receptive field
    spsignals = spcomp' * reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 
    % convolve with temporal filter
    generators(:, iblock) = conv(spsignals, flip(kt), 'valid');
    %======================================================================
    % Diff of Gaussians
    % reduce with spatial receptive field
    spsignalcent = spcompcent' * reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 
    spsignalsurr = spcompsurr' * reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 
    % convolve with temporal filter
    generators(:, iblock) = ...
        conv(spsignalcent, flip(kt), 'valid') + conv(spsignalsurr, flip(ktsurr), 'valid');
    %======================================================================
    % subgrid
    spsignalscent = modelparams.subfilterscent' *...
        reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 
    spsignalssurr =  modelparams.subfilterssurr' *...
        reshape(blockstimulus, [numel(rx) * numel(ry), Nframes]); 

    prenlncent    = conv2(flip(modelparams.kt),     1, spsignalscent', 'valid');
    prenlnsurr    = conv2(flip(modelparams.ktsurr), 1, spsignalssurr', 'valid');

    gens = reshape(logisticfun(prenlncent + prenlnsurr +  modelparams.subbase),...
        size(prenlncent));
    % convolve with temporal filter
    generators(:, iblock,imodel) = gather(gens *  modelparams.subwts);
    %======================================================================
   
end
    
   
generators = reshape(generators, [(Nframes - Nt + 1) * Nblocks, 1]);
%--------------------------------------------------------------------------




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

rxuse = (rx >= xmin & rx <=xmax);
ryuse = (ry >= ymin & ry <=ymax);

xOr = trX - (Nxs / 2) + rx(rxuse);
yOr = trY - (Nys / 2) + ry(ryuse);

xmin = find(rxuse, 1, 'first');
xmax = find(rxuse, 1, 'last');
ymin = find(ryuse, 1, 'first');
ymax = find(ryuse, 1, 'last');

end

