function gmodel = initGaussModelFlicker(sta, cellspikes, stiminfo, orderfit, ...
    spx, spy, ktbas, stascale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
stixelsForFit = 10 * stascale;
rfac   = 10 * 1.4826;

[ypix, xpix, Nt] = size(sta);
[~, imax]    = max(abs(sta(:)));
[y0, x0, t0] = ind2sub(size(sta), imax);
    
rangeX = x0+(-stixelsForFit:stixelsForFit); rangeX = rangeX(rangeX>0 & rangeX<=xpix);
rangeY = y0+(-stixelsForFit:stixelsForFit); rangeY = rangeY(rangeY>0 & rangeY<=ypix);
zoomsta = reshape(sta(rangeY, rangeX, :), numel(rangeY)*numel(rangeX), Nt);

% find significant pixels in the zoomed region
[bpx, ~] = find(abs(zoomsta) > rfac*mad(sta(:),1));

% extract temporal and spatial components
tempcomp = mean(zoomsta(bpx,:),1)';
spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));

baswts = tempcomp' * ktbas;

cfit = abs(spcomp);
cfit = cfit/max(cfit(:));
   
fitprms      = fitgaussrf(spx(rangeX), spy(rangeY), double(cfit));
fitprms(3:4) = fitprms(3:4)/2;

% get rmax
rmax = min(min(fitprms(1), xpix*stascale - fitprms(1)),...
    min(fitprms(2), ypix*stascale - fitprms(2)));

guessactiv = calcGaussianActivationsGrating(fitprms, stiminfo);

allspactivs = guessactiv(orderfit);
allactivs   =  baswts * (ktbas' * allspactivs);
[values, centers] = getNonlinearity(allactivs', cellspikes, 40, 1);
outguess   = fitRLogistic3ToSpikes(double(centers), double(values));
%==========================================================================
% populate model values
gmodel             = struct();
gmodel.ktbasis     = ktbas;
gmodel.ktwts       = baswts * outguess(2);
gmodel.gaussparams = fitprms(1:5);
gmodel.outparams   = [outguess(1) outguess(3)];
gmodel.rmax        = rmax;
%==========================================================================
end

