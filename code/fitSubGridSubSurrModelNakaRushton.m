function [mdlparams, xerr, errall] = fitSubGridSubSurrModelNakaRushton(mdlparams, stiminfo, stimorder, yy, opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
if nargin < 5
    batchsize = 1000;
    eta = 0.01;
    beta1 = .95;
    beta2 = .999;
    epsmall = 1e-8;
    lambda = 1e-3;
    showfig = true;
    Nepochs = 400;
else
    batchsize = opts.batchsize;
    eta = opts.eta;
    beta1 = opts.beta1;
    beta2 = opts.beta2;
    epsmall = opts.epsmall;
    lambda = opts.lambda;
    showfig = opts.showfig;
    Nepochs = opts.Nepochs;
end
%----------------------------------------------------------------------
if ~isfield(mdlparams, 'subwts')
   mdlparams= initializeSubGrid(mdlparams, stiminfo, stimorder, yy);
end
%----------------------------------------------------------------------
if showfig
    fn = figure('Position',[100 100 1700 400]);
    subplot(1,4,1);cla;
    ax = gca; ax.Colormap = flipud(gray); axis equal;
    subplot(1,4,2); cla;
    axis square; ylim([-1 1])
    subplot(1,4,3); cla;
    axis square; xlim([-1 1]); ylim([0 1]);
    subplot(1,4,4); cla;
    axis square; xlim([0 Nepochs])
    xlabel('Epochs'); ylabel('Neg. log-likelihood')
end

xnln    = linspace(-1,1);
Nwt        = numel(mdlparams.ktwts);
Nt         = size(mdlparams.ktbasis, 1);
Nstimuli   = size(stimorder, 2);
Nsubs      = numel(mdlparams.subwts);
pxsize  = mdlparams.pxsize;

iorder = zeros(Nstimuli, Nepochs);
for ii = 1:Nepochs
    iorder(:, ii) = randperm(Nstimuli);
end
iorder = iorder(:);


Nbatches = floor(numel(iorder)/batchsize);

mwold    = gpuArray.zeros(Nsubs + 2*Nwt + 6, 1, 'single');
uwold    = gpuArray.zeros(Nsubs + 2*Nwt + 6, 1, 'single');
xstim    = gpuArray(single(stiminfo));
yspikes  = gpuArray(single(yy));
Np       = sum(yspikes);

errall = gpuArray(NaN(Nbatches, 1, 'single'));
xerr   = (1:Nbatches)*batchsize/Nstimuli;
ArgSubs = cos(calcSubunitActivationPhases(mdlparams, xstim));
finparams = [mdlparams.subwts; mdlparams.ktwts; mdlparams.subsurrwt;...
        mdlparams.subsigma; mdlparams.subsurrsc;...
        mdlparams.subbase; mdlparams.outparams];
    
dij = sqrt((mdlparams.subcnts(:,1) -  mdlparams.subcnts(:,1)').^2 + ...
    (mdlparams.subcnts(:,2) -  mdlparams.subcnts(:,2)').^2);
wd = (dij/min(dij(dij>0))).^-2;
wd(eye(size(wd),'logical')) = 0;
%==========================================================================
% setup bounds
lb  = [zeros(Nsubs,1); -Inf(2 * Nwt,1); 0.7*7.5/pxsize;   1;  -Inf; 0;   0;   0];
ub  = [  Inf(Nsubs,1);  Inf(2 * Nwt,1);   5*7.5/pxsize; Inf;  Inf; Inf; Inf; Inf];
lb  = gpuArray(single(lb));
ub  = gpuArray(single(ub));
%==========================================================================

for ibatch = 1:Nbatches
    
    istart = (ibatch - 1) * batchsize + 1;
    iend   = min(ibatch * batchsize, numel(iorder));
    iuse   = iorder(istart:iend);
    %etause = eta * min(1, 3*ibatch/Nbatches);
    etause = eta * gausslinefun([Nbatches/2 Nbatches/5 1], ibatch);

    [fall, gbatch,nrange] = mleGrad(mdlparams, xstim, ...
        yspikes(iuse), ArgSubs, stimorder(:, iuse), lambda, Np, wd);
    errall(ibatch) = neglogliperspike(yspikes(iuse), fall);
    
    oldparams = [mdlparams.subwts; mdlparams.ktwts; mdlparams.subsurrwt;...
        mdlparams.subsigma;mdlparams.subsurrsc; ...
        mdlparams.subbase; mdlparams.outparams];

    % project gradient
    gproj = boundedProjectedGradient(gbatch, oldparams, lb, ub);
%     gproj = [gbatch(1:Nsubs+2*Nwt);...
%         boundedProjectedGradient(gbatch(Nsubs+2*Nwt+1:end), oldparams(Nsubs+2*Nwt+1:end), lb, ub)];
    % do ADAM
    
    mw = beta1 * mwold +(1-beta1) * gproj;
    uw = beta2 * uwold +(1-beta2) * gproj.^2;
    
    mwold = mw;
    uwold = uw;
    
    mup = mw/(1 - beta1^ibatch);
    uup = uw/(1 - beta2^ibatch);
    
    newparams = oldparams - etause * mup./(sqrt(uup) + epsmall);
    
    % project result
%     newparams = [newparams(1:Nsubs+2*Nwt);...
%         boundedProjectedValues(newparams(Nsubs+2*Nwt+1:end), lb, ub)];
    newparams = boundedProjectedValues(newparams, lb, ub);
        
    finparams = beta1 * finparams +(1-beta1) * newparams;
    
    mdlparams.subwts    = newparams(1:Nsubs);
    %mdlparams.kappas    = newparams(1:Nsubs);
    %mdlparams.subwts    = exp(mdlparams.kappas);
    mdlparams.ktwts     = newparams(Nsubs + (1:Nwt));
    mdlparams.subsurrwt = newparams(Nsubs + Nwt + (1:Nwt));
    mdlparams.subsigma  = newparams(Nsubs + 2 * Nwt + 1);
    mdlparams.subsurrsc = newparams(Nsubs + 2 * Nwt + 2);
    mdlparams.subbase   = newparams(Nsubs + 2 * Nwt + 3);
    mdlparams.outparams = newparams(Nsubs + 2 * Nwt + (4:6));
    %--------------------------------------------------------------------

    %--------------------------------------------------------------------
    if showfig && mod(ibatch, 100) == 1
               
        figure(fn);
        subplot(1,4,1);cla;
        plotSubunitGrid(double(gather(mdlparams.subcnts)),...
            double(gather(mdlparams.subwts)))
        
        title(sprintf('subdiam: %2.1f, sub surr scale: %2.1f', ...
            mdlparams.subsigma*pxsize*4, mdlparams.subsurrsc))
        xlim(mdlparams.gridcent(1) + [-1 1]*50*7.5/pxsize)
        ylim(mdlparams.gridcent(2) + [-1 1]*50*7.5/pxsize)
        
        
        subplot(1,4,2);cla;
        kt     = mdlparams.ktbasis * mdlparams.ktwts;
        surrkt = mdlparams.ktbasis * mdlparams.subsurrwt;
        line(1:Nt, kt/max(abs(kt)), 'Color', 'b')
        line(1:Nt, surrkt/max(abs(kt)), 'Color', 'r')
        
        subplot(1,4,3);cla;
        subparams = [mdlparams.subbase 1];
        xnln      = linspace(min(nrange), max(nrange));
        subv      = rlogistic2(subparams, xnln);
        line(xnln, subv, 'Color', 'k')
        xlim(gather(nrange));
        ylim(gather([min(subv) max(subv)]))

        subplot(1,4,4);cla; 
        errplot = movmean(errall, 100);
        line(xerr, errplot, 'LineStyle', '-');
        
        drawnow;
    end
    %--------------------------------------------------------------------
end

errall = gather(errall);


alldst = sqrt(sum((mdlparams.subcnts-mdlparams.gridcent).^2,2));
iedge  = alldst>max(alldst) - 16/pxsize;
mdlparams.subwts(iedge) = 0;
isilence1 = mdlparams.subwts < max(mdlparams.subwts(~iedge))*.05;
mdlparams.subwts(isilence1)   = 0;
[gparams, ~] = fitCenterGaussSurr(mdlparams);
mdlparams.gaussparams = gparams;
c = getEllipseFromNewParams(mdlparams.gaussparams, 2.5);
isilence2 = ~inpolygon(mdlparams.subcnts(:,1), mdlparams.subcnts(:,2), ...
    c(1,:), c(2,:));
mdlparams.subwts (isilence2)   = 0;
%mdlparams.subwts = mdlparams.subwts/sum(mdlparams.subwts);
isilence = isilence2 | isilence1 | iedge;

[Rfc, Rfs] = calcSubunitSurrActivationAmplitude(mdlparams, stiminfo(:, 1));
kt     = gather(mdlparams.ktbasis * mdlparams.ktwts);
ktsurr = gather(mdlparams.ktbasis * mdlparams.subsurrwt);

Nsubsuse = nnz(~isilence);

subcentacts  = Rfc .* ArgSubs(:, ~isilence);
subcentacts = gather(subcentacts);
tempcentacts  = reshape(subcentacts(stimorder, :), [Nt, Nstimuli * Nsubsuse]);
tempcent      =  kt' * tempcentacts; clear tempcentacts;

subsurracts  = Rfs .* ArgSubs(:, ~isilence);
subsurracts = gather(subsurracts);
tempsurracts = reshape(subsurracts(stimorder, :), [Nt, Nstimuli * Nsubsuse]);
tempsurr     =  ktsurr' * tempsurracts; clear tempsurracts;

allActs  = tempcent + tempsurr;
allActs  = reshape(allActs, [Nstimuli, Nsubsuse]);

allActs  = reshape(logisticfun(allActs + mdlparams.subbase), [Nstimuli, Nsubsuse]);
gensignal = allActs * mdlparams.subwts(~isilence);

outguess            = fitOutputNakaRushton(double(gather(gensignal)), double(yy));
mdlparams.subwts    = mdlparams.subwts;
mdlparams.outparams = gpuArray(single(outguess(:)));

mdlparams.ktbasis   = gpuArray(single(mdlparams.ktbasis));
mdlparams.ktwts     = gpuArray(single(mdlparams.ktwts));
mdlparams.subsurrwt = gpuArray(single(mdlparams.subsurrwt));

if showfig
    figure(fn);
    subplot(1,4,1);cla;
    plotSubunitGrid(double(gather(mdlparams.subcnts)),...
        double(gather(mdlparams.subwts)))
end
%==========================================================================
end
%==========================================================================

function [mdlparams] = initializeSubGrid(mdlparams,xvals, xorder,yy)


dgrid = 16;
pxsize = mdlparams.pxsize;
mdlparams.subsigma = 50/pxsize/4;
mdlparams.cwindow   = 450/pxsize;
mdlparams.subbase   = 0;

% initialize surround diameter
centdiam = getRFDiam(getGaussFromNewParams(mdlparams.gaussparams),2,1);
subdiam = 4*mdlparams.subsigma;
mdlparams.subsurrsc = min(centdiam*mdlparams.surrsc/subdiam,10); 
mdlparams.subsurrwt	= mdlparams.surrktwts;
%----------------------------------------------------------------------
% setup a grid of points that won't change
Nsubs = 1200;
mainpts  = generateHexSubunitGrid(Nsubs);
mdlparams.gridcent = mdlparams.gaussparams(1:2);
gridpts =   mdlparams.gridcent + mainpts * dgrid/pxsize;
mdlparams.subcnts = single(gridpts);
%----------------------------------------------------------------------
% initialize weights and constants properly
subwts = getGaussianSubWeights(mdlparams);
c = getEllipseFromNewParams(mdlparams.gaussparams,2);
subwts(~inpolygon(mdlparams.subcnts(:,1), mdlparams.subcnts(:,2),c(1,:), c(2,:))) = 0;
subwts = subwts/sum(subwts);

[Rfc, Rfs] = calcSubunitSurrActivationAmplitude(mdlparams, xvals(:, 1));
ArgSubs    = cos(calcSubunitActivationPhases(mdlparams, xvals));

Nwt        = numel(mdlparams.ktwts);
Nt         = size(mdlparams.ktbasis, 1);
Nstimuli   = size(xorder, 2);
Nsubsuse   = 50;
[~, isort] = sort(subwts, 'descend');
isubuse    = isort(1:Nsubsuse);

kt    = mdlparams.ktbasis * mdlparams.ktwts';
ktsurr = mdlparams.ktbasis * mdlparams.subsurrwt';

% ktfac = 4 / sum(abs(kt));
% kt = ktfac * kt; 


subcentacts   = Rfc .* ArgSubs(:, isubuse);
tempcentacts  = reshape(subcentacts(xorder, :), [Nt, Nstimuli * Nsubsuse]);
tempcent      =  kt' * tempcentacts; clear tempcentacts;

subsurracts  = Rfs .* ArgSubs(:, isubuse);
tempsurracts = reshape(subsurracts(xorder, :), [Nt, Nstimuli * Nsubsuse]);
tempsurr     =  ktsurr' * tempsurracts; clear tempsurracts;

allActs  = tempcent + tempsurr;
allActs  = reshape(allActs, [Nstimuli, Nsubsuse]);
allActs  = reshape(logisticfun(allActs + mdlparams.subbase), [Nstimuli, Nsubsuse]);

gensignal = allActs * subwts(isubuse);
outguess  = fitOutputNakaRushton(double(gather(gensignal)), double(yy));


mdlparams.subwts    = gpuArray(single(subwts(:)));
mdlparams.outparams = gpuArray(single(outguess(:)));
mdlparams.ktbasis   = gpuArray(single(mdlparams.ktbasis));
mdlparams.ktwts     = gpuArray(single(mdlparams.ktwts))';
mdlparams.subsurrwt = gpuArray(single(mdlparams.subsurrwt))';


end

function [fall, gbatch, nrange] = mleGrad(mdlparams, xvals, yy, ArgSubs, stimorder,lambda, Np, wd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Nt, Nwt] = size(mdlparams.ktbasis);
Nsubs     = numel(mdlparams.subwts);
Nstim     =  size(stimorder,2);
kt        =  mdlparams.ktbasis * mdlparams.ktwts;
ktsurr    = mdlparams.ktbasis * mdlparams.subsurrwt;
    
% forward and backward pass
[Rfc, Rfs, Rjc, Rjs] = calcSubunitSurrActivationAmplitude(mdlparams, xvals(:, 1));
subcentacts  = Rfc .* ArgSubs;
subsurracts  = Rfs .* ArgSubs;

tempcentacts  = reshape(subcentacts(stimorder, :), [Nt, Nstim * Nsubs]);
tempsurracts  = reshape(subsurracts(stimorder, :), [Nt, Nstim * Nsubs]);

tempfin  = kt' * tempcentacts + ktsurr' * tempsurracts;
tempfin  = reshape(tempfin, [size(stimorder,2), Nsubs]);

RsubsNlin   = reshape(logisticfun(tempfin + mdlparams.subbase),[Nstim, Nsubs]);
insignal    = RsubsNlin * mdlparams.subwts;
outparams   = mdlparams.outparams;
[fall,jout] = nakarushton(outparams, insignal);
    
nrange = [min(tempfin,[],'all') max(tempfin,[],'all')];

if nargout > 1
    %----------------------------------------------------------------------
    % triage weights for faster gradient calculations
    cwts     = mdlparams.subwts;
    sumall   = sum(cwts);
    wtvals   = sort(cwts, 'descend');
    iuse     = cwts > wtvals(100) & cwts > 0;
    wtfac    = sumall/sum(cwts(iuse));
    Nsubsuse = nnz(iuse);

    tempactsred      = reshape(subcentacts(stimorder, iuse), [Nt, Nstim * Nsubsuse]);
    tempsurractsred  = reshape(subsurracts(stimorder, iuse), [Nt, Nstim * Nsubsuse]);
    tempbaseprojcent = mdlparams.ktbasis' * tempactsred;
    tempbaseprojsurr = mdlparams.ktbasis' * tempsurractsred;
    %----------------------------------------------------------------------
    % calcuate gradients
    dfdz  = outparams(2) * fall .*(1-fall/outparams(1)) ./insignal;

    dNdg  = reshape(RsubsNlin .* (1 - RsubsNlin), [1, Nstim, Nsubs]);
    dNdgred = dNdg(:, :, iuse);
    
    subjc  = Rjc .* ArgSubs(:, iuse);
    jcacts = reshape(subjc(stimorder, :), [Nt, Nstim * Nsubsuse]);
    jcacts  = reshape(kt' * jcacts, [size(stimorder,2), Nsubsuse]);
    
    subjs  = Rjs .* reshape(ArgSubs(:, iuse), [size(ArgSubs,1), 1, Nsubsuse]);
    jsacts = reshape(subjs(stimorder, :), [Nt, Nstim * 2* Nsubsuse]);
    jsacts  = reshape(ktsurr' * jsacts, [size(stimorder,2),  2, Nsubsuse]);
    
    jsacts(:, 1, :) = squeeze(jsacts(:, 1, :))  + jcacts;
    Jsigma          = (permute(dNdgred, [2 1 3]) .* jsacts) ;
    Jsigma = reshape(Jsigma, [Nstim *2, Nsubsuse]) * mdlparams.subwts(iuse) * wtfac;
    Jsigma = reshape(Jsigma, [Nstim, 2]);
    
    Jbase  = squeeze(dNdg) * mdlparams.subwts;
    
    Jkt = reshape(tempbaseprojcent, [Nwt, Nstim, Nsubsuse]) .* dNdgred;
    Jkt = reshape(Jkt, [Nwt*Nstim, Nsubsuse]) * mdlparams.subwts(iuse) * wtfac;
    Jkt = reshape(Jkt, [Nwt, Nstim])';
    
    Jsurrkt = reshape(tempbaseprojsurr, [Nwt, Nstim, Nsubsuse]) .* dNdgred;
    Jsurrkt = reshape(Jsurrkt, [Nwt*Nstim, Nsubsuse]) * mdlparams.subwts(iuse) * wtfac;
    Jsurrkt = reshape(Jsurrkt, [Nwt, Nstim])';
    %----------------------------------------------------------------------
    % collect gradients
    gall  = [dfdz.* [RsubsNlin Jkt Jsurrkt Jsigma Jbase], jout];
    
    gbatch =  -((yy./fall)' * gall - sum(gall, 1))/Np;
    wuse  = wd*mdlparams.subwts;
    gbatch(1:Nsubs) = gbatch(1:Nsubs) +  2* lambda * wuse'.* ((mdlparams.subwts>0))';
    gbatch = gbatch';
    %----------------------------------------------------------------------
end
%==========================================================================
end
