
% demo_flashed_gratings

% add repository path to MATLAB path
addpath(genpath('..\subunit_grid_model'))
%%
% set experiment folder
expfolder = '20201112_252MEA_mouse_left_half_dorsal';
% load data from selected experiment
expdata      = load(fullfile('..\subunit_grid_model', expfolder, 'expdata.mat'));
% grating flash data
gflashdata   = load(fullfile('..\subunit_grid_model', expfolder, 'gratingflashes_data.mat'));
% imagesequence data
imgseqdata   = load(fullfile('..\subunit_grid_model', expfolder, 'imagesequence_data.mat'));

% set up useful variables
pxsize  = expdata.projector.pixelsize * 1e6; % pixel size in um
screenx = expdata.projector.screen(2); %
screeny = expdata.projector.screen(1);
spX     = 0.5 + (0:screenx-1); 
spY     = 0.5 + (0:screeny-1);

% meta parameters for fitting subunit grid models
% gridspacing = 15/pxsize;
cellWindow  = 800; %in um
% NsubMax     = 2000;
% xvals = linspace(-1, 1);

barwidths = pxsize* sort(0.5./unique(gflashdata.stiminfo(:,1)));
Noris     = numel(unique(gflashdata.stiminfo(:,2)));
Nphases   = numel(unique(gflashdata.stiminfo(:,3)));

pxWindow  = round((cellWindow)/pxsize/2);
xstim     = gflashdata.stiminfo(gflashdata.presentOrder, :);
stimmat   = getGratingMatFromInfo(gflashdata.stiminfo, spX, spY);
Ngratings = size(stimmat, 3);
sizesub   = linspace(-200, 200, 500);


% set up images to match the screen and tranform values into Weber contrast

%imuse = true(size(imgdata.isartificial)); 
imuse          = ~imgseqdata.isartificial; 
imEnsemble     = imgseqdata.imageEnsemble(:, :, imuse);
[Nnaty, Nnatx] = size(imEnsemble, [1 2]);
Nnatimages     = nnz(imuse);
[X, Y] = meshgrid(single(1:screenx), single(1:screeny));
screenImEnsemble = zeros( screeny, screenx, Nnatimages, 'single');
screenImEnsemble((screeny-Nnaty)/2+1:(screeny+Nnaty)/2, ...
    (screenx-Nnatx)/2+1:(screenx+Nnatx)/2, :) = 2*single(imEnsemble)/255-1;
screenImEnsemble = flip(screenImEnsemble, 1);
%%
%--------------------------------------------------------------------------
% select a recorded cell and print its type
icell = 41; icelltypeid = expdata.cellclus_id(icell); %38, 34
if icelltypeid >0,
    typestr = expdata.typelabels{icelltypeid};
    else, typestr = 'Unclassified'; 
end
fprintf('Cell %d selected, type is %s\n', icell, typestr)

cellresp    = gflashdata.trialCounts(icell, :)';
meanresp    = accumarray(gflashdata.presentOrder, cellresp, [Ngratings 1], @mean);
sumresp     = accumarray(gflashdata.presentOrder, cellresp, [Ngratings 1], @sum);
cellrespimg = imgseqdata.trialCounts(icell, :)';
meanrespimg = accumarray(imgseqdata.presentOrder, cellrespimg, [numel(imuse)  1], @mean);

% tuningdata = squeeze(mean(reshape(meanresp, [Nphases, Noris, numel(barwidths)]), 1));
tuningdata = squeeze(max(reshape(meanresp, [Nphases, Noris, numel(barwidths)]), [], 1));

%--------------------------------------------------------------------------
% we will estimate quality measures for our cell, low quality cells may
% show bad model fits/meaningless results
% Symmetrized R2 is not very useful for low-trial data (<3 trials)

% quality for gratings
Ntrialspergrating = accumarray(gflashdata.presentOrder, 1, [Ngratings 1], @sum);
gftrialrsq        = sufTrialRsq( cellresp', gflashdata.presentOrder);
fprintf('Grating symmetrized R2 = %2.2f, mean trials/grating = %2.1f\n', gftrialrsq, mean(Ntrialspergrating))

% quality for images
Ntrialsperimage = accumarray(imgseqdata.presentOrder, 1, [numel(imuse) 1], @sum);
imnattrialrsq   = sufTrialRsq( cellrespimg', imgseqdata.presentOrder);
fprintf('Image   symmetrized R2 = %2.2f, mean trials/image = %2.1f\n', imnattrialrsq, mean(Ntrialsperimage))
%--------------------------------------------------------------------------
% let's first get a gross estimate of the cell's RF to locate its center
csta     = reshape(2*single(stimmat)/255 - 1, screeny*screenx, Ngratings) * single(sumresp);
csta     = reshape(csta, [screeny, screenx])/sum(sumresp);
csta     = csta - mean(csta, 'all');
blursta  = imgaussfilt(csta, 4);
[~, im]  = max(abs(blursta(:)));
[y0, x0] = ind2sub([screeny, screenx], im);

rrx = x0 + (-pxWindow:pxWindow); 
rrx = rrx(rrx > 0 & rrx <= screenx);
rry = y0 + (-pxWindow:pxWindow);
rry = rry(rry > 0 & rry <= screeny);

% center sta and rescale to ease gaussian fit
cfit = csta(rry, rrx);
cfit = abs(cfit);
cfit = cfit - min(cfit,[],'all');
cfit = cfit/max(cfit(:),[],'all');

% perform gaussian fit and calculate 2-sigma ellipse
fitprms = fitgaussrf(spX(rrx), spY(rry), double(cfit));
cel     = getEllipseFromNewParams(fitprms, 2);
%----------------------------------------------------------------------
% We will fit progressively more complex models that are initialized based
% on the previous model in the progression
% MODEL 1: GAUSSIAN RF                + LOGISTIC OUTPUT
% MODEL 2: DIFFERENCE-OF-GAUSSIANS RF + LOGISTIC OUTPUT
% MODEL 3: NONLINEAR SUBUNIT GRID     + NAKA-RUSHTON OUTPUT
%==========================================================================
% MODEL 1: GAUSSIAN RF                + LOGISTIC OUTPUT
%==========================================================================
% initialize
gmdlparams  = struct();
gmdlparams.gaussparams = fitprms(1:5);
rmax = min(min(fitprms(1), screenx - fitprms(1)),...
    min(fitprms(2), screeny - fitprms(2)));

%----------------------------------------------------------------------
% construct initial guess
guessactiv = calcGaussianActivationsGrating(gmdlparams.gaussparams, gflashdata.stiminfo);
outguess   = fitRLogistic3ToSpikes(guessactiv, meanresp);
gmdlparams.outparams = outguess;
gmdlparams.rmax = rmax;

% fit
mdlparams1 = gfFitGaussMonotonic(gmdlparams, xstim, cellresp);

% predict training data
cel1      = getEllipseFromNewParams(mdlparams1.gaussparams, 2);
lfpred    = funfitGaussianModelParams(mdlparams1, gflashdata.stiminfo);    
lfpredll1 = funfitGaussianModelParams(mdlparams1, xstim);    

tuning1 = squeeze(max(reshape(lfpred, [Nphases, Noris, numel(barwidths)]), [], 1));
% tuning1 = squeeze(mean(reshape(lfpred, [Nphases, Noris, numel(barwidths)]), 1));

activ1 = calcGaussianActivationsGrating(mdlparams1.gaussparams, gflashdata.stiminfo);

rsq1 = rsquare(meanresp, lfpred);
nll1 = neglogliperspike(cellresp, lfpredll1);

% predict natural image activations
[imacts1,sta1] = predictGaussModel(mdlparams1, screenImEnsemble, rrx, rry);
rho1           = corr(imacts1', meanrespimg(imuse), 'Type', 'Spearman');

%==========================================================================
% MODEL 2: DIFFERENCE-OF-GAUSSIANS RF + LOGISTIC OUTPUT
%==========================================================================
% fit circ gauss + surround
mdlparams2 = mdlparams1;
mdlparams2.surrsc = 2;
mdlparams2.surrwt = 0.1;

mdlparams2 = gfFitDoGMonotonic(mdlparams2, xstim, cellresp);

cel2   = getEllipseFromNewParams(mdlparams2.gaussparams, 2);
lfpred2   = funfitDoGModelParams(mdlparams2, gflashdata.stiminfo);   
lfpredll2 = funfitDoGModelParams(mdlparams2, xstim);    
tuning2 = squeeze(max(reshape(lfpred2, [Nphases, Noris, numel(barwidths)]), [], 1));
% tuning2 = squeeze(mean(reshape(lfpred2, [Nphases, Noris, numel(barwidths)]), 1));

activ2 = calcDoGActivationsGrating(mdlparams2, gflashdata.stiminfo);

% predict natural image activations

[imacts2, sta2] = predictDogModel(mdlparams2, screenImEnsemble, rrx, rry);
rho2            = corr(imacts2',meanrespimg(imuse), 'Type', 'Spearman');
%==========================================================================
% MODEL 3: NONLINEAR SUBUNIT GRID     + NAKA-RUSHTON OUTPUT
%==========================================================================
% add pixel size in our starting guess
mdlparams2.pxsize = pxsize;

% set optimization parameters
opts         = getDefaultSGparams('flashes');
opts.Nepochs = round(4e5/numel( gflashdata.presentOrder));
opts.showfig = true;

% Obtain a good model initialization
mdlparams3init = gfFitSubunitGridModel(mdlparams2, xstim, cellresp, opts,true);

% Run optimization with a particular regularization strength (lambda)
opts.lambda                = 1e-4;
[mdlparams3, xerr, errfit] = gfFitSubunitGridModel(mdlparams3init, xstim, cellresp, opts,false);


Rsubs          = calcSubunitGridOutputs(mdlparams3, gflashdata.stiminfo);
subunitoutputs = rlogistic2(mdlparams3.subparams, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams3.subwts;
lfpred3        = nakarushton( mdlparams3.outparams, insig);

tuning3 = squeeze(max(reshape(gather(lfpred3), [Nphases, Noris, numel(barwidths)]), [], 1));
% tuning3 = squeeze(mean(reshape(gather(lfpred3), [Nphases, Noris, numel(barwidths)]), 1));

activ3 = subunitoutputs * mdlparams3.subwts/sum(mdlparams3.subwts);
activ3 = gather(activ3);
%outuse = [mdlparams3.bconst sum(mdlparams3.subwts) mdlparams3.outparams];


Rsubs          = calcSubunitGridOutputs(mdlparams3, xstim);
subunitoutputs = rlogistic2(mdlparams3.subparams, Rsubs);

subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams3.subwts;
lfpredll3      = nakarushton( mdlparams3.outparams, insig);

[imacts3, sta3] = predictSingleSubunitModel(mdlparams3, screenImEnsemble, rrx, rry);

rho3  = corr(imacts3',meanrespimg(imuse), 'Type', 'Spearman');
%% summary figure with model fits

% You can navigate all cells and experiments to get an understanding of their
% nonlinear RF properties!
modelstr = {'Gaussian', 'Diff. of Gaussians', 'Subunit grid'};

fw   = 18; fh = 15; % figure width and height in cm
fsub = figure('Color','w','Units', 'centimeters');
fsub.ToolBar  = 'none'; %f1.MenuBar='none';
fsub.Position = [2 2 fw fh]; 
fsub.Renderer = 'painters';

p = panel();
p.pack('v', 4);
for ii = 1:4
    p(ii).pack('h', {0.3 0.3 0.2 0.2})
end
p.fontsize  = 8;
p.de.margin = 1;
p.margin = [1 10 1 10];
p.de.margintop = 10;
for ii = 1:4
   p(ii, 3).marginleft = 6; 
end

p(1,2).select(); cla;
plotLogTuning2D(barwidths, 1:Noris, tuningdata);
caxis([min(tuningdata(:)) max(tuningdata(:))])
yticks([1 Noris]); xticks([15 960])
ylabel('Orientation')

p(1,1).select(); cla;
imagesc(spX(rrx), spY(rry), cfit, [-1 1] * max(abs(cfit), [], 'all'))
ax = gca; ax.Colormap = flipud(gray);
axis equal; xlim(spX([min(rrx) max(rrx)])); ylim(spY([min(rry) max(rry)]))
ax.Visible = 'off'; ax.Title.Visible = 'on';
title('Init. filter (STA of gratings)')

p.title(sprintf('Cell %d, %s\n', icell, typestr))

gractsall  = [activ1 activ2 activ3];
stamodels  = cat(3, sta1, sta2, sta3);
imactsall  = [imacts1' imacts2' imacts3'];
modeltun   = cat(3, tuning1, tuning2, tuning3);
meanimresp = meanrespimg(imuse);
maxyim     = ceil(max(meanimresp)/2) * 2;
maxygr     = ceil(max(meanresp)/2) * 2;

for imodel = 1:3
    
    stacurr = stamodels(:,:,imodel);
    p(1+imodel,1).select(); cla;
    imagesc(spX(rrx), spY(rry), stacurr, [-1 1] * max(abs(stacurr), [], 'all'))
    ax = gca; ax.Colormap = flipud(gray); ax.Visible = 'off';
    ax.Title.Visible = 'on';
    title(modelstr{imodel})
    axis equal; xlim(spX([min(rrx) max(rrx)])); ylim(spY([min(rry) max(rry)]))

    p(1+imodel,2).select();cla;
    plotLogTuning2D(barwidths, 1:Noris, modeltun(:, :, imodel));
    caxis([min(tuningdata(:)) max(tuningdata(:))])
    yticks([1 Noris]); xticks([15 960])
    ylabel('Orientation')
    if imodel==3
        xlabel('Bar width (um)')
    end
    
    
    p(1+imodel,3).select();cla;
    axis square; ylim([0 maxygr]);
    yticks([0 maxygr/2 maxygr])
    line(gractsall(:, imodel), meanresp,...
        'MarkerSize',2,'Marker','o','LineStyle','none')
    ylabel('Spike count')
    title('Grating flashes')
    if imodel==3
        xlabel('Rec. field prediction')
  
    end
    
    p(1+imodel,4).select();cla;
    axis square; ylim([0 maxyim]);
    yticks([0 maxyim/2 maxyim])
    line(imactsall(:, imodel), meanrespimg(imuse),...
        'MarkerSize',3,'Marker','o','LineStyle','none')
    rho = corr(imactsall(:, imodel),meanrespimg(imuse), 'Type', 'Spearman');
    title(sprintf("Nat. images (rho = %2.3f)", rho))
    if imodel==3
        xlabel('Rec. field prediction')
    end
end
%% you can use the code below to save generated figure for the cell
filename = fullfile('..\subunit_grid_model', sprintf('%s_cell_%d.pdf', expfolder, icell));
p.export(filename, sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp')


