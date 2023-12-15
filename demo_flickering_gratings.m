

% demo_flickering_gratings

% add repository path to MATLAB path
addpath(genpath('..\subunit_grid_model'))

%%
% set experiment folder
expfolder = '20210803_252MEA_mouse_left_half_dorsal';
% load data from selected experiment
expdata      = load(fullfile('..\subunit_grid_model', expfolder, 'expdata.mat'));
% grating flash data
gflickerdata = load(fullfile('..\subunit_grid_model', expfolder, 'gratingflicker_data.mat'));
% imagesequence data
fixmovdata   = load(fullfile('..\subunit_grid_model', expfolder, 'fixationmovie_data.mat'));
%%
% let's initialize useful variables
screenfs = expdata.projector.refreshrate;
spX     = gflickerdata.spX; 
spY     = gflickerdata.spY; 

stimpara = gflickerdata.rawdata(1).stimPara;
xmar     = stimpara.lmargin;
ymar     = stimpara.bmargin;

xpix = numel(spX); ypix = numel(spY);

stiminfo  = gflickerdata.stiminfo;
Nt        = size(gflickerdata.ktbas, 1);
orderfit  = gflickerdata.orderfit;
spikesfit = gflickerdata.spikesfit;

frozenorder = orderGratingFlicker(size(stiminfo,1), stimpara.secondseed, stimpara.FrozenFrames);
hfrmat      = uint16(hankel(1:Nt, Nt:stimpara.FrozenFrames));
frozenorder = frozenorder(hfrmat);
%%
icell = 10; icelltypeid = expdata.cellclus_id(icell); %38
if icelltypeid >0
    typestr = expdata.typelabels{icelltypeid};
    else, typestr = 'Unclassified'; 
end
fprintf('Cell %d selected, type is %s\n', icell, typestr)
%%
%-------------------------------------------------------------------------
% estimate cell quality by symmetrized Rsq in the frozen grating part
% first find frozen spikes
Nparts        = numel(gflickerdata.rawdata);
frozenspikes  = cell(Nparts, 1);
runningspikes = cell(Nparts, 1);

for ipart = 1:Nparts
    partdata  = gflickerdata.rawdata(ipart);
    Nframes   = size(partdata.spikesbin, 2);
    spikesbin = partdata.spikesbin(icell, :);
    runningFrames = partdata.stimPara.RunningFrames;
    frozenFrames  = partdata.stimPara.FrozenFrames;
    trialFrames   = runningFrames+frozenFrames;
    
    Ntrials     = floor(Nframes/trialFrames);
    totalFrames = Ntrials*trialFrames;
    
    totalbin   = reshape(spikesbin(1:totalFrames), trialFrames, Ntrials);
    frozenspikes{ipart}  = totalbin(runningFrames+Nt:end,:);
    runningspikes{ipart}  = totalbin(Nt:runningFrames,:);
end
frozenbin   = cat(2, frozenspikes{:});
runningbin  = cat(2, runningspikes{:});
frozenRates = mean(frozenbin, 2)*screenfs;
gftrialrsq  = imageTrialRsq( reshape(frozenbin', [1, size(frozenbin,2), size(frozenbin,1)]));

fprintf('Grating symmetrized R2 = %2.2f\n', gftrialrsq)
%%
%-------------------------------------------------------------------------
% let's first calculate an sta from the gratings
% the code is parallelizing different cells
% this part is typically long and is using the GPU, but the sta can be only
% calculated once

stimmat  = getGratingMatFromInfo(stiminfo, spX, spY);

% we scale down our stimulus for calculating the STA, as it is only used
% for initializing model fitting
scfac = 5;
stimdown = zeros(numel(spY)/scfac,numel(spX)/scfac, size(stimmat,3), 'single');
for ii = 1:size(stimmat,3)
    stimdown(:,:,ii) = imresize(stimmat(:,:,ii),1/scfac,'Method','box','Antialiasing',false);
end

sta = calculateGratingFlickerSTA(reshape(stimdown, [ypix*xpix/scfac^2, size(stiminfo,1)]),...
    gflickerdata.spikesfit(icell, :), gflickerdata.orderfit);
sta = reshape(sta, [ypix/scfac, xpix/scfac, Nt]);

%-------------------------------------------------------------------------
%%
%-------------------------------------------------------------------------
cellspikes    = spikesfit(icell, :)';

cellspikes_nat = squeeze(fixmovdata.runningbin(icell, Nt:end, :));
cellspikes_nat = cellspikes_nat(:);
frmovierates   = squeeze(mean(fixmovdata.frozenbin(icell, Nt:end, :), 3))'*screenfs;
%-------------------------------------------------------------------------
%%
% initialize Gauss model

dspx = mean(reshape(spX, scfac, []), 1);
dspy = mean(reshape(spY, scfac, []), 1);
ingaussmodel = initGaussModelFlicker(sta, cellspikes, stiminfo, orderfit, ...
    dspx, dspy, gflickerdata.ktbas, scfac);
ingaussmodel.screenx = expdata.projector.screen(2);
ingaussmodel.screeny = expdata.projector.screen(1);
ingaussmodel.pxsize = expdata.projector.pixelsize*1e6;


%-------------------------------------------------------------------------
% fit actual Gauss model
fprintf('Fitting Gaussian RF model... ');tic;
mdlparams1 = fitGaussModel(ingaussmodel, stiminfo, ...
    gflickerdata.stimorder, double(cellspikes));
fprintf('Done! Time elapsed: %2.2f s\n', toc);


% get prediction for frozen part of grating movie
pred1 = predictGaussModelFlicker(mdlparams1, stiminfo, frozenorder);
rsq1  = rsquare( frozenRates, pred1*screenfs);
fprintf('Gauss model performance for frozen gratings, R^2 = %2.2f\n', rsq1);

% get prediction for natural movie
natgens1 = generatorsFixationMovie(mdlparams1, ...
    fixmovdata.runningImages, fixmovdata.runningfixations,xmar, ymar);
frozengens1 = generatorsFixationMovie(mdlparams1, ...
    fixmovdata.frozenImages, fixmovdata.frozenfixations, xmar, ymar);

frpreds1 = fitAndPredictOutputNL(natgens1, frozengens1, cellspikes_nat, 'logistic4');
frpreds1 = frpreds1 * screenfs;
natrsq1 = rsquare( frmovierates, frpreds1);
fprintf('Gauss model performance for natural movie, R^2 = %2.2f\n\n', natrsq1);
%-------------------------------------------------------------------------
%%
% fit DoG model
fprintf('Fitting Difference of Gaussians (DoG) RF model... ');tic;
mdlparams2 = fitDoGModel(mdlparams1, stiminfo, gflickerdata.stimorder, double(cellspikes));
fprintf('Done! Time elapsed: %2.2f s\n', toc);

% get prediction for frozen part of grating movie

pred2 = predictDoGModelFlicker(mdlparams2, stiminfo, frozenorder);
rsq2  = rsquare( frozenRates, pred2*screenfs);
fprintf('Diff. of Gaussians model performance for frozen gratings, R^2 = %2.2f\n', rsq2);

% get prediction for natural movie

natgens2 = generatorsFixationMovie(mdlparams2, ...
    fixmovdata.runningImages, fixmovdata.runningfixations,xmar, ymar);
frozengens2 = generatorsFixationMovie(mdlparams2, ...
    fixmovdata.frozenImages, fixmovdata.frozenfixations, xmar, ymar);
        
frpreds2 = fitAndPredictOutputNL(natgens2, frozengens2, cellspikes_nat, 'logistic4');
frpreds2 = frpreds2 * screenfs;
natrsq2 = rsquare( frmovierates, frpreds2);
fprintf('Diff. of Gaussians model performance for natural movie, R^2 = %2.2f\n\n', natrsq2);

%%
% fit Subunit Grid model

opts = struct();

opts.batchsize = 2000;
opts.eta = 0.02;
opts.beta1 = .9;
opts.beta2 = .999;
opts.epsmall = 1e-6;
opts.lambda = 50e-5;
opts.showfig = true;
opts.Nepochs = 40;

fprintf('Fitting Subunit Grid (SG) RF model... ');tic;
mdlparams3 = fitSubGridSubSurrModelNakaRushton(...
    mdlparams2, stiminfo, orderfit, cellspikes, opts);
fprintf('Done! Time elapsed: %2.2f s\n', toc);

% get prediction for frozen part of grating movie

pred3 = predictSubGridSubSurrModel(mdlparams3, stiminfo, frozenorder);
rsq3  = rsquare( frozenRates, pred3*screenfs);
fprintf('Diff. of Gaussians model performance for frozen gratings, R^2 = %2.2f\n', rsq2);

% get prediction for natural movie
natgens3    = generatorsFixationMovie(mdlparams3, ...
    fixmovdata.runningImages, fixmovdata.runningfixations, xmar, ymar);
frozengens3 = generatorsFixationMovie(mdlparams3, ...
    fixmovdata.frozenImages,   fixmovdata.frozenfixations, xmar, ymar);

frpreds3 = fitAndPredictOutputNL(natgens3, frozengens3, cellspikes_nat, 'logistic4');
frpreds3 = frpreds3 * screenfs;
natrsq3 = rsquare( frmovierates, frpreds3);
fprintf('Subunit grid model performance for natural movie, R^2 = %2.2f\n\n', natrsq3);
%--------------------------------------------------------------------------

%%



figure;
subplot(2,1,1)
plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred2*screenfs);
title(sprintf('Frozen grating prediction\n qualityRsq =  %2.2f, RsqDoG = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq2));
ylabel('Firing rate (sp/s)')

subplot(2,1,2)
plot(fixmovdata.frozentimes, fixmovdata.frozenRates(icell, :), ...
    fixmovdata.frozentimes(Nt:end), natpred2)
title(sprintf('Natural movie prediction\n RsqDoG = %2.2f, RsqWN = %2.2f',...
    natrsq2, fixmovdata.fmValidRsq(icell)));
xlim([0 max(fixmovdata.frozentimes)])
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
% scontinds = imds(1:Nstim);
% extrainds = 1:numel(effscnot);
% extrainds(scontinds) =[];
% scont(Nt+imds(1:Nstim))
%%
figure;
subplot(2,1,1)

plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred4*screenfs);
title(sprintf('Frozen grating prediciton\nqualityRsq  =  %2.2f, RsqSubGridSubSurr = %2.3f', ...
    gfdata.allReliableRsq(icell), rsq4));
ylabel('Firing rate (sp/s)')
xlim([0 max(ftimevec)])

subplot(2,1,2)
plot(fixmovdata.frozentimes, fixmovdata.frozenRates(icell, :), ...
    fixmovdata.frozentimes(Nt:end), natpred4)
title(sprintf('Natural movie prediction\nRsqSubGridSubSurr = %2.2f, RsqWN = %2.2f',...
    natrsq4, fixmovdata.fmValidRsq(icell)));
ylabel('Firing rate (sp/s)')
xlabel('Time (s)')
xlim([0 max(fixmovdata.frozentimes)])

figure;
subplot(2, 1, 1)
plot(ftimevec(Nt:end), frozenRates(icell, :), ...
    ftimevec(Nt:end), pred1*screenfs);
title(sprintf('Frozen grating prediction\n qualityRsq =  %2.2f, Rsq = %2.3f',...
    gfdata.allReliableRsq(icell), rsq1));

ylabel('Firing rate (sp/s)')

subplot(2, 1, 2)
plot(fixmovdata.frozentimes, fixmovdata.frozenRates(icell, :), ...
    fixmovdata.frozentimes(Nt:end), natpred1)
title(sprintf('Natural movie prediction\n RsqGauss = %2.2f, RsqWN = %2.2f',...
    natrsq1, fixmovdata.fmValidRsq(icell)));
xlim([0 max(fixmovdata.frozentimes)])
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

