% add repository path to MATLAB path
addpath(genpath('..\subunit_grid_model'))
%%

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

pxWindow = round((cellWindow)/pxsize/2);
xstim    = gflashdata.stiminfo(gflashdata.presentOrder, :);
stimmat  = getGratingMatFromInfo(gflashdata.stiminfo, spX, spY);
Ngratings = size(stimmat, 3);
sizesub  = linspace(-200, 200, 500);


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

% select a recorded cell and print its type
icell = 38; typestr = expdata.typelabels{expdata.cellclus_id(icell)};  
fprintf('Cell %d selected, type is %s\n', icell, typestr)

cellresp   = gflashdata.trialCounts(icell, :)';
meanresp   = accumarray(gflashdata.presentOrder, cellresp, [Ngratings 1], @mean);
%natres     = imgseqdata.meanCounts(icell, imuse);

% let's first get a gross estimate of the cell's RF

csta     = reshape(2*single(stimmat)/255 - 1, screeny*screenx, Ngratings) * single(meanresp);
csta     = reshape(csta, [screeny, screenx])/sum(meanresp);
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
%----------------------------------------------------------------------
% MODEL 1: GAUSSIAN RF                + LOGISTIC OUTPUT
% MODEL 2: DIFFERENCE-OF-GAUSSIANS RF + LOGISTIC OUTPUT
% MODEL 3: NONLINEAR SUBUNIT GRID     + NAKA-RUSHTON OUTPUT

gmdlparams  = struct();
gmdlparams.gaussparams = fitprms(1:5);

% get rmax
rmax = min(min(fitprms(1), screenx - fitprms(1)),...
    min(fitprms(2), screeny - fitprms(2)));

%----------------------------------------------------------------------
% construct initial guess
guessactiv = gfmodels.calcGaussianActivationsGrating(gmdlparams.gaussparams, gfdata.stiminfo);
outguess   = gfmodels.fitRLogistic3ToSpikes(guessactiv, meanresp);
gmdlparams.outparams = outguess;
gmdlparams.rmax = rmax;

p(1,1).select(); cla;
ax = gca; axis equal; ax.Colormap = gray;
xlim(fitprms(1) + [-1 1] * 60); ylim(fitprms(2) + [-1 1] * 60);
imagesc(spX, spY, csta, [-1 1] * max(abs(csta(:))));
line(contpts(1,:), contpts(2,:), 'Color', 'g')
title('Data');
tuningdata = squeeze(mean(reshape(meanresp, [Nphases, Noris, numel(barwidths)]), 1));
%tuningdata = squeeze(max(reshape(respfit, [4, 12, 25]), [],1));

p(1,4).select(); cla;
plotLogTuning2D(barwidths, 1:Noris, tuningdata);
caxis([0 max(tuningdata(:))])

p(1,5).select(); cla;
line(dogacts, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(dogacts', natres', 'Type', 'Spearman'));
title(tstr)
%----------------------------------------------------------------------
% fit gauss + monotonic nonlin
mdlparams1 = gfmodels.gfFitGaussMonotonic(gmdlparams, xstim, cellresp);

cel1      = getEllipseFromNewParams(mdlparams1.gaussparams, 2);
lfpred    = gfmodels.funfitGaussianModelParams(mdlparams1, gfdata.stiminfo);    
lfpredll1 = gfmodels.funfitGaussianModelParams(mdlparams1, xstim);    

tuning1 = squeeze(mean(reshape(lfpred, [Nphases, Noris, numel(barwidths)]), 1));
%tuning1 = squeeze(max(reshape(lfpred, [4, 12, 25]), [], 1));

p(2,1).select(); cla;
ax = gca; axis equal; ax.Colormap = gray;
xlim(fitprms(1) + [-1 1] * 60); ylim(fitprms(2) + [-1 1] * 60);
imagesc(spX, spY, csta, [-1 1] * max(abs(csta(:))));
line(cel1(1,:), cel1(2,:), 'Color', 'r')
line(noiserf(1,:), noiserf(2,:), 'Color', 'g')

title('Gauss Monotonic');

guessactiv = gfmodels.calcGaussianActivationsGrating(mdlparams1.gaussparams, gfdata.stiminfo);

p(2,3).select(); cla;
xlim([-1 1]);
line(guessactiv, meanresp,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
line(xvals, rlogistic3(mdlparams1.outparams,xvals), 'Color','r','Linewidth', 1.5)

p(2,4).select();
plotLogTuning2D(barwidths, 1:Noris, tuning1);
caxis([0 max(tuningdata(:))])
rsq1 = rsquare(meanresp, lfpred);
nll1 = neglogliperspike(cellresp, lfpredll1);
title(sprintf('Rsq = %2.3f, nLL/spike = %2.3f', rsq1, nll1))

imacts1 = predictGaussModel(mdlparams1, screenImEnsemble, rangeX, screeny - rangeY);
imacts1r = imacts1 .*( imacts1 >0);

p(2,5).select(); cla;
line(imacts1, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(imacts1', natres', 'Type', 'Spearman'));
title(tstr)

%----------------------------------------------------------------------
% fit circ gauss + surround
mdlparams2 = mdlparams1;
mdlparams2.surrsc = 2;
mdlparams2.surrwt = 0.1;

mdlparams2 = gfmodels.gfFitDoGMonotonic(mdlparams2, xstim, cellresp);

cel2   = getEllipseFromNewParams(mdlparams2.gaussparams, 2);
lfpred2   = gfmodels.funfitDoGModelParams(mdlparams2, gfdata.stiminfo);   
lfpredll2 = gfmodels.funfitDoGModelParams(mdlparams2, xstim);    
tuning2 = squeeze(mean(reshape(lfpred2, [Nphases, Noris, numel(barwidths)]), 1));

%tuning2 = squeeze(max(reshape(lfpred2, [4, 12, 25]),[], 1));

guessactiv = gfmodels.calcDoGActivationsGrating(mdlparams2, gfdata.stiminfo);

p(3,1).select(); cla;
ax = gca; axis equal; ax.Colormap = gray;
xlim(fitprms(1) + [-1 1] * 60); ylim(fitprms(2) + [-1 1] * 60);
imagesc(spX, spY, csta, [-1 1] * max(abs(csta(:))));
line(cel2(1,:), cel2(2,:), 'Color', 'r')
line(noiserf(1,:), noiserf(2,:), 'Color', 'g')
title('DoG Monotonic');


p(3,3).select(); cla;
xlim([-1 1] * max(abs(guessactiv)));
line(guessactiv, meanresp,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
line(xvals, rlogistic3(mdlparams2.outparams,xvals), 'Color','r','Linewidth', 1.5)


p(3,4).select();
plotLogTuning2D(barwidths, 1:Noris, tuning2);
caxis([0 max(tuningdata(:))])
rsq2 = rsquare(meanresp, lfpred2);
nll2 = neglogliperspike(cellresp, lfpredll2);
title(sprintf('Rsq = %2.3f, nLL/spike = %2.3f', rsq2, nll2))

imacts2 = predictDogModel(mdlparams2, screenImEnsemble);

p(3,5).select(); cla;
xlim([-1 1] * max(abs(imacts2)));
line(imacts2, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(imacts2', natres', 'Type', 'Spearman'));
title(tstr)
mdlparams2.pxsize = pxsize;

%----------------------------------------------------------------------
%%
ops.batchsize = 64;
ops.eta = 0.01;
ops.beta1 = .9;
ops.beta2 = .999;
ops.epsmall = 1e-6;
ops.lambda = 1e-4;
ops.showfig = true;
ops.Nepochs = round(4e5/size(gfdata.trialCounts,2));
ops.delaylambda = false;
ops.delaysubunit = false;
% 
% tic;
% [mdlparams3, ~, ~] = gfmodels.gfFitRank(mdlparams2, gfdata.stiminfo, respfit, ops);
% toc;

% tic;
% [mdlparams3, ~, ~] = gfmodels.gfFitNonlinearSubunitGapL1(...
%     mdlparams2, xstim, cellresp, ops);
% toc;
tic;
[mdlparams3, ~, ~] = gfmodels.gfFitNonlinearSubunitGapL1swish(...
    mdlparams2, xstim, cellresp, ops,false);
toc;
%%

ops.lambda = 1e-2;
tic;
mdlparams3 = gfmodels.gfFitSubGridBasis3(mdlparams2, xstim, cellresp,ops);
toc;

% mdlarray = gather([mdlparams3.subwts; mdlparams3.subsigma; mdlparams3.subsurrsc;...
%     mdlparams3.subsurrwt; mdlparams3.subparams; mdlparams3.bconst; mdlparams3.outparams]);
% calcCoverageFactor(mdlarray',mdlparams3.gridcent)
% tic;
% [mdlparams3, ~, ~] = gfmodels.gfFitNonlinearSubunitGapMonoOut(mdlparams2, gfdata.stiminfo, respfit, ops);
% toc;

% tic;
% %mdlparams3 = gfmodels.gfFitNonlinearSubunitSingleMleFast(mdlparams2, xstim, cellresp);
% mdlparams3 = gfmodels.gfFitNonlinearSubunitSingle2(mdlparams2, xstim, cellresp);
% toc;
%%
Rsubs          = gfmodels.calcSubunitGridOutputsExp(mdlparams3, gfdata.stiminfo);
%subunitoutputs = rlogistic2(mdlparams3.subparams, Rsubs);
subunitoutputs = subnonlinbasis(mdlparams3, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams3.subwts + mdlparams3.bconst;
% lfpred3        = mdlparams3.outparams./(1 + exp(-insig));
lfpred3        = nakarushton( mdlparams3.outparams, insig);

%tuning3 = squeeze(max(reshape(lfpred3, [4, 12, 25]),[], 1));
tuning3 = squeeze(mean(reshape(lfpred3, [Nphases, Noris, numel(barwidths)]), 1));

activ3 = subunitoutputs * mdlparams3.subwts/sum(mdlparams3.subwts);
activ3 = gather(activ3);
%outuse = [mdlparams3.bconst sum(mdlparams3.subwts) mdlparams3.outparams];


Rsubs          = gfmodels.calcSubunitGridOutputsExp(mdlparams3, xstim);
%subunitoutputs = rlogistic2(mdlparams3.subparams, Rsubs);
subunitoutputs = subnonlinbasis(mdlparams3, Rsubs);

subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams3.subwts + mdlparams3.bconst;
%lfpredll3      = mdlparams3.outparams./(1 + exp(-insig));
lfpredll3      = nakarushton( mdlparams3.outparams, insig);

[imrf,pxranges] = generateSubunitImage(mdlparams3, fitprms(1:2), 60);

p(4,1).select(); cla;
axis square;
xlim(fitprms(1) + [-1 1] * 60); ylim(fitprms(2) + [-1 1] * 60);
imagesc(pxranges(:,1), pxranges(:,2), imrf,[0 1])
line(contpts(1,:), contpts(2,:), 'Color', 'g')
ax =gca; ax.Colormap = flipud(gray);
title('Subunit Grid Monotonic')
% 
p(4,2,1).select(); cla;
axis square; xlim([-1 1] * 130); ylim([-0.2 1])
outplot =  gfmodels.linedogplot(mdlparams3.subsigma * pxsize,...
        mdlparams3.subsurrsc, mdlparams3.subsurrwt, sizesub);
line(sizesub, outplot); xlim([min(sizesub) max(sizesub)]);
xlabel('X position (um)')
title(sprintf('Sub. diam = %d um', round(4 * mdlparams3.subsigma * pxsize)))

p(4,2,2).select(); cla;
axis square;
nlnplot  = rlogistic2(mdlparams3.subparams, xvals);
nlnplot  = nlnplot/max(nlnplot);
line(xvals, nlnplot); xlim([-1 1]); ylim([0 1]);
xlabel('Subunit input'); ylabel('Subunit output')

% 
xout = linspace(min(activ3), max(activ3));
p(4,3).select(); cla;
xlim([min(activ3) max(activ3)]);
line(activ3, meanresp,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
%line(xout, rlogistic3(outuse,xout), 'Color','r','Linewidth', 1.5)
line(xout, nakarushton(mdlparams3.outparams,xout), 'Color','r','Linewidth', 1.5)

% 
p(4,4).select(); cla;
plotLogTuning2D(barwidths, 1:Noris, tuning3);
caxis([0 max(tuningdata(:))])
rsq3 = rsquare(meanresp, lfpred3);
nll3 = neglogliperspike(cellresp, lfpredll3);
title(sprintf('Rsq = %2.3f, nLL/spike = %2.3f', rsq3, nll3))
%%
imacts3 = predictSingleSubunitModel(mdlparams3, screenImEnsemble, rangeX, screeny - rangeY);

p(4,5).select(); cla;
xlim(gather([min(imacts3) max(imacts3)]));
line(imacts3, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(imacts3', natres', 'Type', 'Spearman'));
title(tstr)
% 
% [~, graySuf ]= sufImg( sufdata.stimPara, screenx, screeny);
% [cta, ctb, Ncombis, prs] = findSubunitFlashOrder(sufdata.stimPara, 1000);
% 
% sufImEnsemble = zeros(screeny, screenx, numel(cta),'single');
% for ii = 1:numel(cta)
%     currim = zeros(screeny, screenx,'single');
%     currim(graySuf<0) = ctb(ii);
%     currim(graySuf>0) = cta(ii);
%     sufImEnsemble(:, :, ii) = currim;
% end
% 
% sufImEnsemble = flip(sufImEnsemble, 1);
% sufacts = predictSingleSubunitModel(mdlparams3, sufImEnsemble, rangeX, screeny - rangeY);
% sufinsig = sufacts*sum( mdlparams3.subwts) + mdlparams3.bconst;
% sufpreds  = mdlparams3.outparams./(1 + exp(-sufinsig));
% xcon = linspace(-1, 1, 2*sufdata.stimPara.ncontrasts +1);
% 
% figure;
% sufresp = reshape(sufpreds,17,17);
% subplot(1,3,1); 
% imagesc(xcon,xcon,sufresp)
% ax = gca; ax.YDir = 'normal'; 
% axis equal; xlim([-1 1]); ylim([-1 1])
% subplot(1,3,2);
% contour(xcon,xcon,sufresp)
% ax = gca; ax.YDir = 'normal'; 
% axis equal; xlim([-1 1]); ylim([-1 1])
% subplot(1,3,3);
% plot(sufinsig, sufdata.meanCounts(icell,:),'o')
% sufcorr = corr(sufinsig', sufdata.meanCounts(icell,:)', 'Type','Spearman');
% title(sprintf('corr: %2.3f', sufcorr))

% savepath = 'C:\Users\Karamanlis_Dimokrati\Documents\DimosFolder\conferences\202107_retinal circuits symposium\poster';
% filename = 'example_onoffcell.pdf';
% p.export(fullfile(savepath,filename), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp')
% % 

%%
ops.lambda  = 1e-4;

tic;%
[mdlparams4,~] = gfmodels.gfFitNonlinearSubunitDoubleGapL1(mdlparams3, xstim, cellresp, ops);
toc;

Rsubson    = gfmodels.calcSubunitGridOutputs(mdlparams4.onparams, gfdata.stiminfo);
suboutson  = rlogistic2(mdlparams4.onparams.subparams, Rsubson);
suboutson  = reshape(suboutson, size(Rsubson));
Rsubsoff   = gfmodels.calcSubunitGridOutputs(mdlparams4.offparams, gfdata.stiminfo);
suboutsoff = rlogistic2(mdlparams4.offparams.subparams, Rsubsoff);
suboutsoff = reshape(suboutsoff, size(Rsubson));
insig      = suboutson * mdlparams4.onparams.subwts + ...
   suboutsoff * mdlparams4.offparams.subwts + mdlparams4.bconst;
lfpred4    = mdlparams4.outparams./(1 + exp(-insig));

tuning4 = squeeze(mean(reshape(lfpred4, [4, Noris, numel(barwidths)]), 1));

sum4   = sum(mdlparams4.onparams.subwts) + sum(mdlparams4.offparams.subwts);
activ4 = (insig - mdlparams4.bconst)/sum4;
activ4 = gather(activ4);
outuse = [mdlparams4.bconst sum4 mdlparams4.outparams];



Rsubson    = gfmodels.calcSubunitGridOutputs(mdlparams4.onparams, xstim);
suboutson  = rlogistic2(mdlparams4.onparams.subparams, Rsubson);
suboutson  = reshape(suboutson, size(Rsubson));
Rsubsoff   = gfmodels.calcSubunitGridOutputs(mdlparams4.offparams, xstim);
suboutsoff = rlogistic2(mdlparams4.offparams.subparams, Rsubsoff);
suboutsoff = reshape(suboutsoff, size(Rsubson));
insig      = suboutson * mdlparams4.onparams.subwts + ...
   suboutsoff * mdlparams4.offparams.subwts + mdlparams4.bconst;
lfpredll4   = mdlparams4.outparams./(1 + exp(-insig));


xout = linspace(min(activ4), max(activ4));
p(5,3).select(); cla;
xlim([min(activ4) max(activ4)]);
line(activ4, meanresp,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
line(xout, rlogistic3(outuse,xout), 'Color','r','Linewidth', 1.5)


p(5,4).select(); cla;
plotLogTuning2D(barwidths, 1:Noris, tuning4);
caxis([0 max(tuningdata(:))])
rsq4 = rsquare(meanresp, lfpred4);
nll4 = neglogliperspike(cellresp, lfpredll4);
title(sprintf('Rsq = %2.3f, nLL/spike = %2.3f', rsq4, nll4))

oox = round((range(rangeX))/2);
ooy = round((range(rangeY))/2);

newRangeX = round(mdlparams4.onparams.gridcent(1)) + (-oox:oox);
newRangeY = round(mdlparams4.onparams.gridcent(2)) + (-ooy:ooy);

imacts4 = predictDoubleSubunitModel(mdlparams4, screenImEnsemble, newRangeX, newRangeY);
imacts4 = gather(imacts4);


p(5,5).select(); cla;
xlim([min(imacts4) max(imacts4)]);
line(imacts4, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(imacts4', natres', 'Type', 'Spearman'));
title(tstr)



%%
tic;
mdlparams4 = gfmodels.gfFitNonlinearSubunitUshapeMleFastMahal(mdlparams3, xstim, cellresp);
%mdlparams3 = gfmodels.gfFitNonlinearSubunitSingleMleGauss(mdlparams2, xstim, cellresp);
toc;

Rsubs          = gfmodels.calcSubunitGridOutputs(mdlparams4, gfdata.stiminfo);
subunitoutputs = ulogistic5gpu(mdlparams4.subparams, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams4.subwts + mdlparams4.bconst;
lfpred4        = mdlparams4.outparams./(1 + exp(-insig));

%tuning3 = squeeze(max(reshape(lfpred3, [4, 12, 25]),[], 1));
tuning4 = squeeze(mean(reshape(lfpred4, [4, Noris, numel(barwidths)]), 1));

activ4 = subunitoutputs * mdlparams4.subwts/sum(mdlparams4.subwts);
activ4 = gather(activ4);
outuse = [mdlparams4.bconst sum(mdlparams4.subwts) mdlparams4.outparams];

Rsubs          = gfmodels.calcSubunitGridOutputs(mdlparams4, xstim);
subunitoutputs = ulogistic5gpu(mdlparams4.subparams, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs));
insig          = subunitoutputs * mdlparams4.subwts + mdlparams4.bconst;
lfpredll4      = mdlparams4.outparams./(1 + exp(-insig));

[imrf,pxranges] = generateSubunitImage(mdlparams4, fitprms(1:2), 60);

p(5,1).select(); cla;
axis square;
xlim(fitprms(1) + [-1 1] * 60); ylim(fitprms(2) + [-1 1] * 60);
imagesc(pxranges(:,1), pxranges(:,2), imrf,[0 1])
line(contpts(1,:), contpts(2,:), 'Color', 'g')
ax =gca; ax.Colormap = flipud(gray);
title('Subunit Grid Monotonic')
% 
p(5,2,1).select(); cla;
axis square;
outplot =  gfmodels.linedogplot(mdlparams4.subsigma * pxsize,...
        mdlparams4.subsurrsc, mdlparams4.subsurrwt, sizesub);
line(sizesub, outplot); xlim([min(sizesub) max(sizesub)]);
xlabel('X position (um)')
title(sprintf('Sub. diam = %d um', round(4 * mdlparams4.subsigma * pxsize)))

p(5,2,2).select(); cla;
axis square;
nlnplot  = ulogistic5(mdlparams4.subparams, xvals);
nlnplot  = nlnplot/max(nlnplot);
line(xvals, nlnplot); xlim([-1 1]); ylim([0 1]);
xlabel('Subunit input'); ylabel('Subunit output')

% 
xout = linspace(min(activ4), max(activ4));
p(5,3).select(); cla;
xlim([min(activ4) max(activ4)]);
line(activ4, meanresp,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
line(xout, rlogistic3(outuse,xout), 'Color','r','Linewidth', 1.5)

% 
p(5,4).select(); cla;
surf(barwidths, 1:Noris, tuning4, 'EdgeColor', 'none'); view(2);
caxis([0 max(tuningdata(:))])
rsq4 = rsquare(meanresp, lfpred4);
nll4 = neglogliperspike(cellresp, lfpredll4);
title(sprintf('Rsq = %2.3f, nLL/spike = %2.3f', rsq4, nll4))


imacts4 = predictUSubunitModel(mdlparams4, screenImEnsemble, rangeX, screeny - rangeY);
imacts4 = gather(imacts4);

p(5,5).select(); cla;
xlim([min(imacts4) max(imacts4)]);
line(imacts4, natres,'Color','k','Marker','o', 'LineStyle','none','MarkerSize',3); 
tstr = sprintf('Spearman rho = %2.3f', corr(imacts4', natres', 'Type', 'Spearman'));
title(tstr)

sufacts = predictUSubunitModel(mdlparams4, sufImEnsemble, rangeX, screeny - rangeY);
sufinsig = sufacts*sum( mdlparams4.subwts) + mdlparams4.bconst;
sufpreds  = mdlparams4.outparams./(1 + exp(-sufinsig));
xcon = linspace(-1, 1, 2*sufdata.stimPara.ncontrasts +1);

figure;
sufresp = reshape(sufpreds,17,17);
subplot(1,3,1); 
imagesc(xcon,xcon,sufresp)
ax = gca; ax.YDir = 'normal'; 
axis equal; xlim([-1 1]); ylim([-1 1])
subplot(1,3,2);
contour(xcon,xcon,sufresp)
ax = gca; ax.YDir = 'normal'; 
axis equal; xlim([-1 1]); ylim([-1 1])
subplot(1,3,3);
plot(sufinsig, sufdata.meanCounts(icell,:),'o')
sufcorr = corr(sufinsig', sufdata.meanCounts(icell,:)', 'Type','Spearman');
title(sprintf('corr: %2.3f', sufcorr))



