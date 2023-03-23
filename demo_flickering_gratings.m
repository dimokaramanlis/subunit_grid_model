%%
%experiment = callExperimentFromList(88, true);
experiment = callExperimentFromList(88, true);

%experiment = callExperimentFromList(72, true);

gfdata     = loadAnalyzedStimData(experiment, 'gratingflicker');
gfastimdata = loadRawStimData(experiment, experiment.stimids.gratingflicker(1));
% imId   = experiment.stimids.imagesequence(1);
% imgdata = loadAnalyzedStimData(experiment, imId, ...
%     {'sprankNatural', 'imageEnsemble', 'isartificial','meanCounts', ...
%     'spaceVecX', 'spaceVecY', 'allActivations'});

%-------------------------------------------------------------------------
% load movie data

fmid = experiment.stimids.fixationmovie(1);
fixdata =  loadAnalyzedStimData(experiment, fmid, ...
    {'options', 'frozenRates', 'frozentimes', 'fmValidRsq','linearpreds'});
stimdata = loadRawStimData(experiment,fmid);
fprintf('Separating running from frozen data... '); tic;
[frozenbin, frozenRasters, runningbin, frozenfixations, runningfixations] = ...
    loadFixationStimuliAndSpikesMouse(experiment, stimdata, fixdata.options);

[idstoload, ~, newids] = unique(runningfixations(1,:));
imEnsemble = loadSingleFrameImEnsemble(stimdata.stimPara, idstoload);

% rewrite runnning list with the new ids
runningfixations(1, :) = newids;

[frozenidstoload, ~, newfrozenids] = unique(frozenfixations(1,:));
imFrozenEnsemble = loadSingleFrameImEnsemble(stimdata.stimPara, frozenidstoload);

% rewrite runnning list with the new ids
frozenfixations(1, :) = newfrozenids;
fprintf('Done! Took %2.2f s\n', toc);
%-------------------------------------------------------------------------

chid = experiment.stimids.frozencheckerflicker(1);
chdata = loadAnalyzedStimData(experiment, chid, ...
    {'contourpoints','rfmodelparams','temporalComponents','timeVec'});

screenfs = experiment.projector.refreshrate;
pxsize  = experiment.projector.pixelsize*1e6;
screenx = experiment.projector.screen(2);
screeny = experiment.projector.screen(1);
spX = gfdata.spX; spY = gfdata.spY; 

xmar = gfastimdata.stimPara.lmargin;
ymar = gfastimdata.stimPara.bmargin;

xpix = numel(spX); ypix = numel(spY);

gridspacing = 16/pxsize;
cellWindow  = 800; %in um
NsubMax     = 2000;
xvals = linspace(-1, 1);

pxWindow = round((cellWindow)/pxsize/2);
mainpts  = generateHexSubunitGrid(NsubMax);

stiminfo  = gfdata.stiminfo;
Nt        = size(gfdata.staAll, 3);
orderfit  = gfdata.orderfit;
spikesfit = gfdata.spikesfit;
frozenorder = gfdata.frozenorder;
frozenRates = gfdata.frozenRates;
ftimevec    = gfdata.frozenTimeVec;

% stnmfpath  = fullfile(experiment.originalFolder,'data_analysis', 'stnmf');
% stnmfdata = load(fullfile(stnmfpath, 'stnmf.mat'));
% stnmfids  = [stnmfdata.cells(:).cluster_id];


%%
icell      = 21; close all;

cellspikes = spikesfit(icell, :)';

natcellspikes = squeeze(runningbin(icell, Nt:end, :));
natcellspikes = natcellspikes(:);
        
%-------------------------------------------------------------------------
% % ----- Represent this filter (approximately) in temporal basis -------
ktbasprs.neye = min(2,Nt); % Number of "identity" basis vectors near time of spike;
ktbasprs.ncos = min(8,Nt); % Number of raised-cosine vectors to use  
ktbasprs.kpeaks = [0 ((Nt-ktbasprs.neye)/2)];  % Position of 1st and last bump 
ktbasprs.b = 1; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(ktbasprs,Nt);
%-------------------------------------------------------------------------

stixelsForFit = 50;
rfac   = 10 * 1.4826;

% extract components for center
csta    = reshape(squeeze(gfdata.staAll(icell, :, :)), [ypix, xpix, Nt]);
csta    = 2 * single(csta)/255 - 1;
% fit center surround rf, get center and surround temporal components

[~, imax]    = max(abs(csta(:)));
[y0, x0, t0] = ind2sub(size(csta), imax);
    
rangeX = x0+(-stixelsForFit:stixelsForFit); rangeX = rangeX(rangeX>0 & rangeX<=xpix);
rangeY = y0+(-stixelsForFit:stixelsForFit); rangeY = rangeY(rangeY>0 & rangeY<=ypix);
zoomsta = reshape(csta(rangeY, rangeX, :), numel(rangeY)*numel(rangeX), Nt);

% find significant pixels in the zoomed region
[bpx, ~] = find(abs(zoomsta) > rfac*mad(csta(:),1));

% extract temporal and spatial components
tempcomp = mean(zoomsta(bpx,:),1)';
spcomp   = reshape(zoomsta*tempcomp, numel(rangeY),numel(rangeX));

baswts = tempcomp' * ktbas;

cfit = abs(spcomp);
cfit = cfit/max(cfit(:));
   
fitprms      = fitgaussrf(spX(rangeX), spY(rangeY), double(cfit));
fitprms(3:4) = fitprms(3:4)/2;

% get rmax
rmax = min(min(fitprms(1), xpix - fitprms(1)),...
    min(fitprms(2), ypix - fitprms(2)));

guessactiv = gfmodels.calcGaussianActivationsGrating(fitprms, stiminfo);

allspactivs = guessactiv(orderfit);
allactivs   =  baswts * (ktbas' * allspactivs);
[values, centers] = getNonlinearity(allactivs', cellspikes, 40, 1);
outguess   = gfmodels.fitRLogistic3ToSpikes(double(centers), double(values));

gmdlparams  = struct();
gmdlparams.ktbasis     = ktbas;
gmdlparams.ktwts       = baswts * outguess(2);
gmdlparams.gaussparams = fitprms(1:5);
gmdlparams.outparams   = [outguess(1) outguess(3)];
gmdlparams.rmax        = rmax;

%iuse = randperm(numel(cellspikes), min(numel(cellspikes)/2, 150000));

tic;
mdlparams1 = gflicker.fitGaussModel(gmdlparams, stiminfo, ...
    gfdata.stimorder, double(cellspikes));
toc;

% tic;
% % mdlparams1 = gflicker.fitGaussModel(gmdlparams, stiminfo, ...
% %     orderfit(:, iuse), double(cellspikes(iuse)));
% mdlparams1s = gflicker.fitGaussModel(gmdlparams, stiminfo, ...
%     orderfit, double(cellspikes));
% toc;


pred1 = gflicker.predictGaussModel(mdlparams1, stiminfo, frozenorder);
rsq1 = rsquare( frozenRates(icell, :)', pred1*screenfs);
mdlparams1.ktwts = mdlparams1.ktwts';

% get all activations from rfs and corresponding spikes
natgens = getGaussGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, mdlparams1,  xmar, ymar);
[svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
guess            = fitRLogisticToSpikes(double(scents), double(svals));
nlnparams        = fitOutputNonlinearityML4(double(natgens), double(natcellspikes),  [0 guess(2:end)]);
frozengens1 = getGaussGeneratorsFixMovie(experiment, ...
            imFrozenEnsemble, frozenfixations, mdlparams1, xmar, ymar);
natpred1      = rlogistic4(nlnparams, frozengens1) * screenfs;
%natpred1      = interp1(scents, svals, frozengens1,'linear','extrap')* screenfs;
natrsq1 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred1);

figure;
subplot(2, 1, 1)
plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred1*screenfs);
title(sprintf('Frozen grating prediction\n qualityRsq =  %2.2f, Rsq = %2.3f',...
    gfdata.allReliableRsq(icell), rsq1));

ylabel('Firing rate (sp/s)')

subplot(2, 1, 2)
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred1)
title(sprintf('Natural movie prediction\n RsqGauss = %2.2f, RsqWN = %2.2f',...
    natrsq1, fixdata.fmValidRsq(icell)));
xlim([0 max(fixdata.frozentimes)])
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')

figure;

subplot(1, 2, 1)
c = getEllipseFromNewParams(mdlparams1.gaussparams, 2);
axis equal; 
xlim(mdlparams1.gaussparams(1) + 50 *[-1 1]* 7.5/pxsize); 
ylim(mdlparams1.gaussparams(2) + 50 *[-1 1]* 7.5/pxsize); 
line(c(1,:), c(2, :), 'Color', 'k')
mdlparams1.ktwts = mdlparams1.ktwts';

subplot(1, 2, 2)
axis square; xlim([1 Nt])
line(1:Nt, mdlparams1.ktwts*mdlparams1.ktbasis', 'Color', 'k')

%%

tic;
% mdlparams2 = gflicker.fitDoGModel(mdlparams1, stiminfo, ...
%     orderfit(:, iuse), double(cellspikes(iuse)));
mdlparams2 = gflicker.fitDoGModel(mdlparams1, stiminfo, ...
    gfdata.stimorder, double(cellspikes));
toc;


pred2      = gflicker.predictDoGModel(mdlparams2, stiminfo, frozenorder);
rsq2 = rsquare( frozenRates(icell, :)', pred2*screenfs);

% % get all activations from rfs and corresponding spikes
mdlparams2.ktwts = mdlparams2.ktwts';
mdlparams2.surrktwts = mdlparams2.surrktwts';

natgens = getDoGGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, mdlparams2,  xmar, ymar);
[svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
guess            = fitRLogisticToSpikes(double(scents), double(svals));
nlnparams        = fitOutputNonlinearityML4(double(natgens), double(natcellspikes),  [0 guess(2:end)]);
frozengens2 = getDoGGeneratorsFixMovie(experiment, ...
            imFrozenEnsemble, frozenfixations, mdlparams2, xmar, ymar);
natpred2      = rlogistic4(nlnparams, frozengens2) * screenfs;
%natpred2      = interp1(scents, svals, frozengens2,'linear','extrap')* screenfs;

natrsq2 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred2);

figure;
subplot(2,1,1)
plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred2*screenfs);
title(sprintf('Frozen grating prediction\n qualityRsq =  %2.2f, RsqDoG = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq2));
ylabel('Firing rate (sp/s)')

subplot(2,1,2)
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred2)
title(sprintf('Natural movie prediction\n RsqDoG = %2.2f, RsqWN = %2.2f',...
    natrsq2, fixdata.fmValidRsq(icell)));
xlim([0 max(fixdata.frozentimes)])
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')

figure;

subplot(1, 2, 1)
c     = getEllipseFromNewParams(mdlparams2.gaussparams, 2);
csurr = getEllipseFromNewParams(mdlparams2.gaussparams, 2 * mdlparams2.surrsc);

axis equal; 
xlim(mdlparams1.gaussparams(1) + 50 *[-1 1] * 7.5/pxsize); 
ylim(mdlparams1.gaussparams(2) + 50 *[-1 1]* 7.5/pxsize); 
line(c(1,:), c(2, :), 'Color', 'b')
line(csurr(1,:), csurr(2, :), 'Color', 'r')
mdlparams2.ktwts = mdlparams2.ktwts';
mdlparams2.surrktwts = mdlparams2.surrktwts';

subplot(1, 2, 2)
axis square; xlim([1 Nt])
line(1:Nt, mdlparams2.ktwts*mdlparams1.ktbasis', 'Color', 'b')
line(1:Nt, mdlparams2.surrktwts*mdlparams1.ktbasis', 'Color', 'r')
mdlparams2.pxsize = pxsize;
%%

opts = struct();

opts.batchsize = 2000;
opts.eta = 0.02;
opts.beta1 = .9;
opts.beta2 = .999;
opts.epsmall = 1e-6;
opts.lambda = 50e-5;
opts.showfig = true;
opts.Nepochs = 40;

% tic;
% mdlparams4 = gflicker.fitSubGridSubSurrModelStart2(...
%     mdlparams2, stiminfo, orderfit, cellspikes, opts);
% toc;

tic;
mdlparams4 = gflicker.fitSubGridSubSurrModelNakaRushton(...
    mdlparams2, stiminfo, orderfit, cellspikes, opts);
toc;

%%
% tic;
% mdlparams3 = gflicker.fitSubGridSubSurrModelBCSC(...
%     mdlparams2, stiminfo, orderfit(:,iuse), cellspikes(iuse), opts);
% toc;

tic;
mdlparams3 = gflicker.fitSubGridModelFull(...
    mdlparams2, stiminfo, gfdata.stimorder, cellspikes, opts);
toc;


pred3 = gflicker.predictSubGridModel(mdlparams3, stiminfo, frozenorder);
rsq3  = rsquare( frozenRates(icell, :)', pred3*screenfs);


%--------------------------------------------------------------------------
% movie preds
% natgens = getSubGridSubSurrGeneratorsFixMovie(experiment, ...
%     imEnsemble, runningfixations, mdlparams3,  xmar, ymar);
natgens = getSubGridGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, mdlparams3,  xmar, ymar);
[svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
guess            = fitRLogisticToSpikes(double(scents), double(svals));
guess(1) = 0;
nlnparams        = fitOutputNonlinearityML4(double(natgens), double(natcellspikes),  guess);
frozengens3 = getSubGridGeneratorsFixMovie(experiment, ...
            imFrozenEnsemble, frozenfixations, mdlparams3, xmar, ymar);
natpred3      = rlogistic4(nlnparams, frozengens3) * screenfs;
%natpred3      = interp1(scents, svals, frozengens3,'linear','extrap')* screenfs;

natrsq3 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred3);
%--------------------------------------------------------------------------


figure;
subplot(2,1,1)

plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred3*screenfs);
title(sprintf('Frozen grating prediciton\nqualityRsq  =  %2.2f, RsqSubGrid = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq3));
ylabel('Firing rate (sp/s)')

subplot(2,1,2)
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred3)
title(sprintf('Natural movie prediction\nRsqSubGrid = %2.2f, RsqWN = %2.2f',...
    natrsq3, fixdata.fmValidRsq(icell)));
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')


%%
opts = struct();

opts.batchsize = 512;
opts.eta = 0.01;
opts.beta1 = .9;
opts.beta2 = .999;
opts.epsmall = 1e-6;
opts.lambda = 1e-6;
opts.showfig = true;
opts.Nepochs = 20;

tic;
mdlparams4 = gflicker.fitSubGridSubSurrModel(...
    mdlparams3, stiminfo, orderfit, cellspikes, opts);
toc;
%%
pred4      = gflicker.predictSubGridSubSurrModel(mdlparams4, stiminfo, frozenorder);
rsq4 = rsquare( frozenRates(icell, :)', pred4*screenfs);


%--------------------------------------------------------------------------
% movie preds
tic;
natgens = getSubGridSubSurrGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, mdlparams4,  xmar, ymar);
toc;

% [svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
% guess            = fitNakaRushtonToSpikes(double(scents), double(svals));
nlnparamsnr        = fitOutputNakaRushton(double(natgens), double(natcellspikes));

% nlnparams        = fitOutputNonlinearityML4(double(natgens), double(natcellspikes),  nguess);
%nlnparams        = fitSoftPlus3ML(double(natgens), double(natcellspikes));

[frozengens4] = getSubGridSubSurrGeneratorsFixMovie(experiment, ...
            imFrozenEnsemble, frozenfixations, mdlparams4, xmar, ymar);
% natpred4      = rlogistic4(nlnparams, frozengens4) * screenfs;
natpred4      = nakarushton(nlnparamsnr, frozengens4) * screenfs;
%natpred4      = softplusfun3m(nlnparams, frozengens4) * screenfs;
%natpred4      = interp1(scents, svals, frozengens4,'linear','extrap')* screenfs;

natrsq4 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred4);
%--------------------------------------------------------------------------
Nstim = round(numel(natpred4)/5);
[aa, imds] = sort(abs(natpred4 - natpred2),'descend');
dataresp =fixdata.frozenRates(icell, Nt:end)';
%rsquare( dataresp(imds(1:Nstim)), natpred2(imds(1:Nstim)))

%effscnot = conv(scont,flip(mdlparams4.ktbasis*mdlparams4.ktwts),'valid');


% scontinds = imds(1:Nstim);
% extrainds = 1:numel(effscnot);
% extrainds(scontinds) =[];
% scont(Nt+imds(1:Nstim))

figure;
subplot(2,1,1)

plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred4*screenfs);
title(sprintf('Frozen grating prediciton\nqualityRsq  =  %2.2f, RsqSubGridSubSurr = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq4));
ylabel('Firing rate (sp/s)')
xlim([0 max(ftimevec)])

subplot(2,1,2)
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred4)
title(sprintf('Natural movie prediction\nRsqSubGridSubSurr = %2.2f, RsqWN = %2.2f',...
    natrsq4, fixdata.fmValidRsq(icell)));
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')
xlim([0 max(fixdata.frozentimes)])
%%
xprof = linspace(-200,200,500);
wtratio = norm(mdlparams4.surrktwts)/norm(mdlparams4.ktwts);
outplot =  gfmodels.linedogplot(mdlparams4.subsigma*2.5,...
    mdlparams4.subsurrsc, wtratio,...
    xprof);

% sigmas = linspace(0,10);
% activations = (1-exp(-sigmas.^2/2))-...
%     wtratio.*(1-exp(-sigmas.^2/2./mdlparams4.subsurrsc.^2));
% [~, im] = max(activations);
% effsigma = sigmas(im);
% 
% plot(linspace(0,500,500),activations)
%%

stid = find(experiment.clusters(icell,1)==stnmfids);
celldata = stnmfdata.cells(stid);
stmods = celldata.modules;
xrange = celldata.crop{1};
xrange = (double(xrange{1}+1:xrange{2})-0.5)*6+0.5;
yrange = celldata.crop{2};
yrange = 600-(double(yrange{1}+1:yrange{2})-0.5)*6+0.5;
[xx, yy] = meshgrid(xrange, yrange);
Nmods  = nnz(celldata.morans_i.subunit>0.25);
gparamsall = NaN(Nmods, 6);
newmods    = zeros(Nmods, numel(yrange), numel(xrange));
for ii = 1:Nmods
    currsub = double(squeeze(stmods(ii,:,:)))';
    [~,im] = max(abs(currsub(:)));
    if currsub(im)< 0
        currsub = -currsub;
    end
    currsub = currsub - min(currsub(:));
    currsub = currsub/max(currsub(:));
    gparamsall(ii,:) = fitgaussrf(xrange,yrange, currsub);
    subnew = gauss2dfun(gparamsall(ii,:), {xx(:), yy(:)});
    subnew = reshape(subnew, numel(yrange), numel(xrange));
    newmods(ii, :, :) = subnew/norm(subnew(:));
%     c = getEllipseFromNewParams(gparamsall(ii,:),2);
%     line(c(1,:), 600-c(2,:))
end
stnmfparams.modules  = newmods;
stnmfparams.tempcomp = chdata.temporalComponents(icell, :);
stnmfparams.Nblinks  = 4;
stnmfparams.Npx      = 6;
stnmfparams.xrange   = xrange;
stnmfparams.yrange   = yrange;
stnmfparams.subwts   = lsqnonneg(double(newmods(:,:))', celldata.sta.space(:));

frefresh = experiment.projector.refreshrate;
mtvec    = (-(Nt-1/2):1:-1/2)*1/frefresh; %in seconds

natgens = getStnmfGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, stnmfparams, chdata.timeVec, mtvec);
[svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
guess            = fitRLogisticToSpikes(double(scents(~isnan(svals))), double(svals(~isnan(svals))));
nguess = [0 guess(2:end)];
nlnparams        = fitOutputNonlinearityML4(double(natgens), double(natcellspikes),  nguess);
frozengens5 = getStnmfGeneratorsFixMovie(experiment, ...
    imFrozenEnsemble, frozenfixations, stnmfparams, chdata.timeVec, mtvec);
natpred5      = rlogistic4(nlnparams, frozengens5) * screenfs;
natrsq5 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred5);

figure;
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred5)
title(sprintf('Natural movie prediction\nRsqSTNMF = %2.2f, RsqWN = %2.2f',...
    natrsq5, fixdata.fmValidRsq(icell)));
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')

%%

figure;
subplot(1,4,1)
plot(xprof, outplot);
axis square; ylim([-0.2 1])
xlabel('X pos. (um)')
title('Subunit RF')

subplot(1,4,2)
plotSubGridCircles(gather(mdlparams4.subcnts), gather(mdlparams4.subwts), gather(mdlparams4.subsigma))
cpts = squeeze(chdata.contourpoints(icell, : ,:));
line(cpts(1,:), 600-cpts(2,:), 'Color','g')
axis equal; axis tight;
xlim(nanmedian(cpts(1,:))+[-50 50]);ylim(600-nanmedian(cpts(2,:))+[-50 50])
title('Subunit grid layout')
%line([410 430], 290*[1 1],'Color','k')
%text(410, 295, '50 um')

subplot(1,4,3)
line(cpts(1,:), 600-cpts(2,:), 'Color','g')
for ii = 1:Nmods
    c = getEllipseFromNewParams(gparamsall(ii,:),2);
    line(c(1,:), 600-c(2,:),'Color','r')
end
axis equal; axis tight;
xlim(nanmedian(cpts(1,:))+[-50 50]);ylim(600-nanmedian(cpts(2,:))+[-50 50])
title('STNMF layout')

subplot(1,4,4)
plotSubGridCircles(gather(mdlparams4.subcnts), gather(mdlparams4.subwts), gather(mdlparams4.subsigma))
for ii = 1:Nmods
    c = getEllipseFromNewParams(gparamsall(ii,:),2);
    line(c(1,:), 600-c(2,:),'Color','r')
end
axis equal; axis tight;
xlim(nanmedian(cpts(1,:))+[-50 50]);ylim(600-nanmedian(cpts(2,:))+[-50 50])
title('overlay')
%%
% ptimes = repmat(fixdata.frozentimes,[2 1])+[-1;1]*mean(diff(fixdata.frozentimes))/2;
% ptimes = [ptimes(1); ptimes(:); ptimes(end)];
% frmean    = fixdata.frozenRates(icell, :);
% 
% figure;
% patch(ptimes,[0;kron(frmean',ones(2,1));0],[1 1 1]*0.4,'EdgeColor','none')
% line(fixdata.frozentimes(Nt:end),fixdata.linearpreds(icell,:),'LineWidth', 1.2,...
%     'Color','g')
% line(fixdata.frozentimes(Nt:end),natpred4, 'Color','r','LineWidth', 1.2)
% 
% 

%%
%%==
opts = struct();

opts.batchsize = 2048;
opts.eta = 0.01;
opts.beta1 = .99;
opts.beta2 = .999;
opts.epsmall = 1e-8;
opts.lambda = 1e-6;
opts.showfig = true;
opts.Nepochs = 30;

mdlsurrguess = mdlparams4;
surrgauss      = mdlparams2.gaussparams;
surrgauss(3:4) = surrgauss(3:4) * mdlparams2.surrsc;
mdlsurrguess.surrgauss = surrgauss';
mdlsurrguess.surrktwts = mdlparams2.surrktwts';

tic;
mdlparams42 = gflicker.fitSubGridWithSurrModel(...
    mdlsurrguess, stiminfo, orderfit, cellspikes, opts);
toc;

pred42      = gflicker.predictSubGridWithSurr(mdlparams42, stiminfo, frozenorder);
rsq42 = rsquare( frozenRates(icell, :)', pred42*screenfs);


%--------------------------------------------------------------------------
% movie preds
natgens = getSubGridWithSurrGeneratorsFixMovie(experiment, ...
    imEnsemble, runningfixations, mdlparams42,  xmar, ymar);
[svals, scents]  = getNonlinearity(natgens, natcellspikes, 40, 1);
guess            = fitRLogisticToSpikes(double(scents), double(svals));
nlnparams        = fitOutputNonlinearityML(double(natgens), double(natcellspikes), guess(2:end));

frozengens42 = getSubGridWithSurrGeneratorsFixMovie(experiment, ...
            imFrozenEnsemble, frozenfixations, mdlparams42, xmar, ymar);
%natpred42      = rlogistic3(nlnparams, frozengens) * screenfs;
natpred42      = interp1(scents, svals, frozengens42)* screenfs;
natrsq42 = rsquare( fixdata.frozenRates(icell, Nt:end)', natpred42);
%--------------------------------------------------------------------------

figure;
subplot(2,1,1)

plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred42*screenfs);
title(sprintf('Frozen grating prediciton\nqualityRsq  =  %2.2f, RsqSubGridGlobalSurr = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq42));
ylabel('Firing rate (sp/s)')

subplot(2,1,2)
plot(fixdata.frozentimes, fixdata.frozenRates(icell, :), ...
    fixdata.frozentimes(Nt:end), natpred42)
title(sprintf('Natural movie prediction\nRsqSubGridGlobalSurr = %2.2f, RsqWN = %2.2f',...
    natrsq42, fixdata.fmValidRsq(icell)));

ylabel('Firing rate (sp/s)')
xlabel('Time (s)')

