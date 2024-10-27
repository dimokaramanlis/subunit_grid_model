function [mdlparams, xerr, errall] = gfFitSubunitGridModel(mdlparams, xvals, yy, opts,initialize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
if nargin < 4
    batchsize = 1000;
    eta = 0.01;
    beta1 = .9;
    beta2 = .999;
    epsmall = 1e-6;
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

if ~isfield(mdlparams, 'subwts')
    mdlparams = initializeSubGrid(mdlparams,xvals, yy);
end

if initialize, return, end

pxsize = mdlparams.pxsize;
Nsubs = numel(mdlparams.subwts);
%----------------------------------------------------------------------
if showfig
    % initialize plots
    fn = figure('Position',[100 100 1700 400]);
    subplot(1,5,1);
    ax = gca; ax.Colormap = flipud(gray); axis equal;
    title('Current subunit grid and weights')
    subplot(1,5,2);
    ax = gca; axis square; 
    title('Fitting progress')
    ylabel('Neg. log-likelihood')
    xlabel('Epochs')
    
    subplot(1,5,3);
    ax = gca; axis square; 
    title('Current subunit nonlinearity')
    ylabel('Output')
    xlabel('Input')
    
    subplot(1,5,4);
    ax = gca; axis square; 
    title('Current subunit profile')
    ylabel('Output')
    xlabel('X position (um)')
    
    subplot(1,5,5);
    ax = gca; axis square; 
    ylabel('Spike count')
    xlabel('Pooled subunit output')
    title('Current output nonlinearity')

end

xnln    = linspace(-1,1);
sizesub = linspace(-150, 150, 1000);

Ndata   = numel(yy);
rng(1);

iorder = zeros(Ndata, Nepochs);
for ii = 1:Nepochs
    iorder(:, ii) = randperm(Ndata);
end
iorder = iorder(:);


Nbatches = floor(numel(iorder)/batchsize);

mwold = gpuArray.zeros(Nsubs + 9, 1, 'single');
uwold = gpuArray.zeros(Nsubs + 9, 1, 'single');
xstim = gpuArray(single(xvals));
yuse  = gpuArray(single(yy));

errall = gpuArray(NaN(Nbatches, 1, 'single'));
xerr   = gpuArray(NaN(Nbatches, 1, 'single'));

ArgSubs = cos(calcSubunitActivationPhases(mdlparams, xvals));

dij = sqrt((mdlparams.subcnts(:,1) -  mdlparams.subcnts(:,1)').^2 + ...
    (mdlparams.subcnts(:,2) -  mdlparams.subcnts(:,2)').^2);
wd = (dij/min(dij(dij>0))).^-2;
wd(eye(size(wd),'logical')) = 0;


alldst = sqrt(sum((mdlparams.subcnts-mdlparams.gridcent).^2,2));
iedge  = alldst>max(alldst) - 16/pxsize;

[~,~,iun] = unique(xvals,'rows');
yavg = accumarray(iun, yy, [], @mean);
%==========================================================================
% setup bounds

lb  = [zeros(Nsubs,1); 0.8*7.5/pxsize;   1; 0; -Inf; -Inf;   0;  0;   0;  0];
ub  = [  Inf(Nsubs,1);   5*7.5/pxsize; Inf; 1;  Inf;  Inf; Inf; Inf; Inf; Inf];
lb  = gpuArray(single(lb));
ub  = gpuArray(single(ub));
%==========================================================================
etainit = eta;

for ibatch = 1:Nbatches
    
    istart = (ibatch - 1) * batchsize + 1;
    iend   = min(ibatch * batchsize, numel(iorder));
    iuse   = iorder(istart:iend);
    
    %etause = eta * min(1, 5*ibatch/Nbatches);
    etause = eta * gausslinefun([Nbatches/2 Nbatches/5 1], ibatch);
    
    
%     if ibatch < ceil(0.2*Nepochs*Ndata/batchsize) && delaylambda
%         luse = 0;
%     else
%         luse = lambda;
%     end
%     
    [~, ~,gbatch] = mleGradSubunitGridModel(mdlparams, ...
        xstim, yuse, ArgSubs, iuse, lambda, wd);
    %errall(ibatch) = neglogliperspike(yuse, fall);
    
    oldparams = [mdlparams.subwts; mdlparams.subsigma; mdlparams.subsurrsc;...
        mdlparams.subsurrwt; mdlparams.subparams'; mdlparams.outparams;mdlparams.outbase];
    
    % project gradient
    gproj = boundedProjectedGradient(gbatch, oldparams, lb, ub);

    % do ADAM
    
    mw = beta1 * mwold +(1-beta1) * gproj;
    uw = beta2 * uwold +(1-beta2) * gproj.^2;
    
    mwold = mw;
    uwold = uw;
    
    mup = mw/(1 - beta1^ibatch);
    uup = uw/(1 - beta2^ibatch);
    
    % parameter update
    newparams = oldparams - etause * mup./(sqrt(uup) + epsmall);
   
    % project result
    newparams = boundedProjectedValues(newparams, lb, ub);

    mdlparams.subwts   = newparams(1:Nsubs);
    mdlparams.subsigma = newparams(Nsubs + 1);
    mdlparams.subsurrsc = newparams(Nsubs + 2);
    mdlparams.subsurrwt = newparams(Nsubs + 3);
    mdlparams.subparams = newparams(Nsubs + (4:5))';
    mdlparams.outparams = newparams(Nsubs + (6:8));
    mdlparams.outbase   = newparams(Nsubs + 9);
    %--------------------------------------------------------------------
    if nargout > 1 && mod(ibatch, 30) == 1
        [fall,insignal]= mleGradSubunitGridModel(mdlparams, xstim, yuse, ArgSubs, 1:numel(yuse), lambda);
        errall(ibatch) = neglogliperspike(yuse, fall);
        xerr(ibatch) = ibatch;
    end
    %--------------------------------------------------------------------
    if showfig && mod(ibatch, 60) == 1
               
        figure(fn);
        subplot(1,5,1);cla;
        plotSubunitGrid(double(gather(mdlparams.subcnts)),...
            double(gather(mdlparams.subwts)))
%         subalphas = double(gather(mdlparams.subwts));
%         plotSubGridCircles(double(gather(mdlparams.subcnts)),...
%             subalphas/max(subalphas), mdlparams.subsigma,'k')
        xlim(mdlparams.gridcent(1) + [-1 1]*50*7.5/pxsize)
        ylim(mdlparams.gridcent(2) + [-1 1]*50*7.5/pxsize)

        subplot(1,5,2);cla; xlim([1 Nepochs])
        line(xerr(~isnan(xerr))*batchsize/Ndata, errall(~isnan(xerr)), 'LineStyle', '-');
        
        subplot(1,5,3);cla;
        subvals = rlogistic2(mdlparams.subparams, xnln);

        subvals = subvals/max(subvals);
        line(xnln, subvals,'Color', 'b');

        subplot(1,5,4);cla;
        outplot =  linedogplot(mdlparams.subsigma*pxsize,...
            mdlparams.subsurrsc, mdlparams.subsurrwt, sizesub);
        line(sizesub, outplot);
        tstr1 = sprintf('Current subunit profile (diam = %2.1f um)',...
            mdlparams.subsigma*7.5*4);
        tstr2 = sprintf('surr. scale = %2.2f, weight = %2.2f',...
            mdlparams.subsurrsc, mdlparams.subsurrwt / mdlparams.subsurrsc^2);
        title({tstr1, tstr2})
        drawnow;

        subplot(1,5,5);cla;
        line(accumarray(iun,gather(insignal),[],@mean), yavg,'Color','k','Marker','o','LineStyle','none');
        line(insignal, nakarushton(mdlparams.outparams,insignal) + mdlparams.outbase,...
            'Color','r','Marker','.','LineStyle','none');
        drawnow;
        
    end

end



%final cleanup
oldsum = sum(mdlparams.subwts);


% plot(mdlparams.subcnts(:,1), mdlparams.subcnts(:,2),'bo', ...
%     mdlparams.subcnts(iedge,1), mdlparams.subcnts(iedge,2),'ro')
% 
% 
% 
mdlparams.subwts(iedge) = 0;

isilence1 = mdlparams.subwts < max(mdlparams.subwts(~iedge))*.05;
mdlparams.subwts(isilence1)   = 0;


[gparams, ~] = fitCenterGaussSurr(mdlparams);
mdlparams.gaussparams = gparams;
c = getEllipseFromNewParams(mdlparams.gaussparams, 2.5);
isilence2 = ~inpolygon(mdlparams.subcnts(:,1), mdlparams.subcnts(:,2), ...
    c(1,:), c(2,:));

mdlparams.subwts (isilence2)   = 0;
mdlparams.subwts = oldsum*mdlparams.subwts/sum(mdlparams.subwts);

isilence = isilence2 | isilence1 | iedge;
% refit output

Rf           = calcSubunitActivationAmplitudes(mdlparams, xvals(:, 1));
subacts      = Rf .* ArgSubs(:, ~isilence);
RsubsNlin    = rlogistic2gpu(mdlparams.subparams, subacts);
RsubsNlin    = reshape(RsubsNlin, size(subacts));
insignal     = RsubsNlin * mdlparams.subwts(~isilence) + mdlparams.bconst;

%guess = [mdlparams.bconst oldsum mdlparams.outparams];
% guess    = double(gather(guess));

insignal = double(gather(insignal));
spikes    = double(yy);


if any(isnan(insignal)) || all(insignal==0)
    return;
end


% fitprms   = gfmodels.fitOutputNLswish(insignal, spikes, guess);
fitprms =  fitOutputNakaRushtonBase(insignal, spikes, ...
    double(gather([mdlparams.outparams;mdlparams.outbase])));

mdlparams.subwts = mdlparams.subwts;
mdlparams.outparams = fitprms(1:3);
mdlparams.outbase   = fitprms(4);

% mdlparams.subwts = mdlparams.subwts * fitprms(2);
% mdlparams.bconst = fitprms(1);
% mdlparams.outparams = fitprms(3);
%==========================================================================

fall           = mleGradSubunitGridModel(mdlparams, xstim, yuse, ArgSubs, 1:numel(yuse), lambda);
errall(ibatch) = neglogliperspike(yuse, fall);
xerr(ibatch) = ibatch;
xerr   = gather(xerr(~isnan(xerr))*batchsize/Ndata);
errall = gather(errall(~isnan(errall)));

%==========================================================================
if showfig
    figure(fn);
    subplot(1,5,1);cla;
    plotSubunitGrid(double(gather(mdlparams.subcnts)),...
        double(gather(mdlparams.subwts)))
end
%==========================================================================
end


function mdlparams = initializeSubGrid(mdlparams,xvals, yy)
% setup grid
pxsize = mdlparams.pxsize;
dgrid = 16;
mdlparams.subsigma  = 50/pxsize/4;
mdlparams.subsurrwt = 1e-1;
mdlparams.surrgauss = mdlparams.gaussparams;
mdlparams.surrwt = 0;
mdlparams.bconst = 0;
% surrdiam = getRFDiam(getGaussFromNewParams(mdlparams.gaussparams),1,1)*mdlparams.surrsc;
% subsurrsc = max(surrdiam/(2*mdlparams.subsigma),2);
% subsurrsc = min(subsurrsc,20);
% mdlparams.subsurrsc = subsurrsc;
mdlparams.subsurrsc = 6;
mdlparams.outbase = 0;
% %----------------------------------------------------------------------
% % assign subunit nonlinearity parameters

if mdlparams.outparams(2) > 0
    mdlparams.subparams  = gpuArray(single([0 4]));
    %mdlparams.subparams  = gpuArray(single([-3 6]));
    %mdlparams.subparams  = gpuArray(single([-1 8]));
else
    %mdlparams.subparams  = gpuArray(single([0 -4]));
    mdlparams.subparams  = gpuArray(single([-3 -6]));
    %mdlparams.subparams  = gpuArray(single([-1 -8]));
end

mdlparams.cwindow   = 450/pxsize;
%----------------------------------------------------------------------
% setup a grid of points that won't change
Nsubs = 1200;
mainpts  = generateHexSubunitGrid(Nsubs);
mdlparams.gridcent = mdlparams.gaussparams(1:2);
gridpts =   mdlparams.gridcent + mainpts * dgrid/pxsize;
mdlparams.subcnts = gpuArray(single(gridpts));
%----------------------------------------------------------------------
sigmas = mdlparams.gaussparams(3:4);
if max(sigmas)/min(sigmas) > 10
    mdlparams.gaussparams(3:4) = mean(sigmas)*0.75/2;
else
    mdlparams.gaussparams(3:4) =sigmas*0.75;
end

subwts = getGaussianSubWeights(mdlparams);
c = getEllipseFromNewParams(mdlparams.gaussparams,2);
subwts(~inpolygon(mdlparams.subcnts(:,1), mdlparams.subcnts(:,2),c(1,:), c(2,:))) = 0;
subwts = subwts/sum(subwts);

Rsubs          = calcSubunitGridOutputs(mdlparams, xvals);
subunitoutputs = rlogistic2(mdlparams.subparams, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs)) * subwts;


fitprms =  fitOutputNakaRushton(double(gather(subunitoutputs)), double(yy));

% if fitprms(2) < 3
%     fitprms(2) = 3;
% end
mdlparams.subwts = subwts;
mdlparams.outparams = gpuArray(single(fitprms))';

% %----------------------------------------------------------------------

% fitprms = fitWeights(modelToRAM(mdlparams),xvals, yy);
% mdlparams.subwts = gpuArray(single(fitprms));

%----------------------------------------------------------------------

% % assign subunit nonlinearity parameters
% 
fitprms = fitNonlinearities(modelToRAM(mdlparams),xvals, yy);
mdlparams.outparams = gpuArray(single(fitprms(3:end-1)));
mdlparams.subparams = gpuArray(single(fitprms(1:2)))';
mdlparams.outbase  = gpuArray(single(fitprms(end)));
% if mdlparams.outparams(2) < 3
%     mdlparams.outparams(2) = 3;
% end



% oldparams = fitprms(1:2);
% xnln = linspace(-1,1);
% newparams = oldparams;
% if max(rlogistic2(oldparams,xnln))<0.5
%     newparams(2) = newparams(2) *2;
% end
% wtfac = max(rlogistic2(oldparams,xnln))/max(rlogistic2(newparams,xnln));
% mdlparams.subparams = newparams';
% mdlparams.subwts = mdlparams.subwts*wtfac;


end


function fitprms = fitNonlinearities(mdlparams,xvals, yy)
mdlfit = mdlparams;

isilence = mdlfit.subwts == 0;

mdlfit.subwts (isilence) = [];
mdlfit.subcnts(isilence,:) = [];

subacts   = calcSubunitGridOutputs(mdlfit, xvals);
%==========================================================================
foptim = @(p) optimizeSubParams(p, subacts, yy, mdlfit);
%==========================================================================
% setup guess from current params
guess  = [mdlfit.subparams'; mdlfit.outparams; mdlfit.outbase];
%guess = [mdlfit.subparams'; 10;2;10];
%==========================================================================
lb    = [-Inf, -10,      0,  0,   0,    0];
ub    = [Inf,    10,   Inf, Inf, Inf, Inf];
%==========================================================================

options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients',false,...
    'HessianFcn', [], 'FiniteDifferenceType', 'central',...
    'MaxIterations',50);
[fitprms,res] = fmincon(foptim, guess, [],[], [], [], lb, ub, [], options);

end


function [f, g, H] = optimizeSubParams(p, xx, yy, otherparams)
params = otherparams;
% params.outparams(4)    = p(end);
% params.subparams = p(1:end-1);
params.subparams = p(1:2);
params.outparams = p(3:end-1);
params.outbase   = p(end);

Np = sum(yy);

switch nargout
    case {0, 1}
        lf = funfitsubparams(params, xx);
    case 2
        [lf, lg] = funfitsubparams(params, xx);
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
    case 3
        [lf, lg] = funfitsubparams(params, xx);
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
        H = ((yy./(lf.^2)).*lg)' * lg/Np;
end
f = -(log(lf)'*yy - sum(lf))/Np; %objective f


end

function [fout, J] = funfitsubparams(params, subacts)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% get subunit centers & weights
subwts    = params.subwts;
outparams = params.outparams;
subparams = params.subparams;
Ns        = numel(subparams);
%==========================================================================
[RsubsNlin, jnlin] = rlogistic2(subparams, subacts);
insignal           = reshape(RsubsNlin, size(subacts)) * subwts;
% f = modifiedswish(outparams, insignal);
%f = shahoutput(outparams, insignal);
[f,jout] = nakarushton(outparams, insignal);
fout = f + params.outbase;
%==========================================================================
if nargout > 1    

%     expuse = exp(-(outparams(3) * insignal + outparams(4)));
%     dfdz  = outparams(2) * f./insignal + outparams(3) * f .*expuse./(1+expuse);

    % shahoutput
    %dfdz  = f .* (outparams(1)./insignal - outparams(2)./(outparams(2)*insignal+1));

    %naka rushton
    dfdz  = outparams(2) * f .*(1-f/outparams(1)) ./insignal;
    
    Jsub  = permute(reshape(jnlin, [size(subacts) Ns]), [1 3 2]);
    Jsub  = reshape(Jsub, size(subacts,1)*Ns, size(subacts,2)) * subwts;
    Jsub  = reshape(Jsub, size(subacts,1), Ns);
    
    %J =  [dfdz.*Jsub f.*expuse./(1+expuse)];
    J =  [dfdz.*Jsub jout ones(size(dfdz))];
end
end


function outweights = fitWeights(mdlparams, xvals, yy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
mdlfit = mdlparams;
Rsubs          = gfmodels.calcSubunitGridOutputs(mdlfit, xvals);

subunitoutputs = rlogistic2(mdlfit.subparams, Rsubs);
subunitoutputs = reshape(subunitoutputs, size(Rsubs));

Nsubs = numel(mdlfit.subwts);
kappas = log(mdlfit.subwts);
kappas(isinf(kappas)) = min(kappas(~isinf(kappas)))*100;
%==========================================================================
% setup guess from current params
%guess = [mdlfit.kappas; mdlfit.outparams(4)];
guess = kappas;
foptim = @(p) optimizeGridWeights(p, subunitoutputs, yy, mdlfit);

options = optimoptions('fminunc','Algorithm','trust-region',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients',false,...
    'HessianFcn','objective','MaxIterations',20);
%==========================================================================
[fitprms, res] = fminunc(foptim, guess, options);

outweights = exp(fitprms);

end





function [f, g, H] = optimizeGridWeights(wtparams, xx, yy, otherparams)

params = otherparams;
% params.subwts = exp(wtparams(1:end-1)); 
% params.kappas = wtparams(1:end-1); 
% params.outparams(4)    = wtparams(end);

params.subwts = exp(wtparams); 
params.kappas = wtparams; 


[lf, lg] = wtToResp(params, xx);
Np       = sum(yy);

f = -(log(lf)'*yy - sum(lf))/Np;%objective f

if nargout > 1
    g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
end

if nargout > 2
    H = ((yy./(lf.^2)).*lg)' * lg/Np;
    %H(1:end-1, 1:end-1) = H(1:end-1, 1:end-1) + lambda * diag(params.subwts);
end
end


function [f, J] = wtToResp(params, xvals)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

outparams = params.outparams;
Aout = params.outparams;
% combine outputs 

fin = xvals * params.subwts;
f = nakarushton(outparams, fin);

if nargout > 1
    %naka rushton
    dfdz  = outparams(2) * f .*(1-f/outparams(1)) ./fin;
    
    J = dfdz .* (xvals .* params.subwts');
end

end
