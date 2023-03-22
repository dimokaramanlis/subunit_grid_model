function params = fitNakaRushtonToSpikes(xvals, spikes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
%unwrap and standardize variables
spikes = spikes(:); 
xx     = xvals(:);
yy     = spikes/max(spikes);
%==========================================================================

guess = [1 1 0.1];

lb = [0 0 0];
ub = [Inf  Inf  Inf];
%==========================================================================
%optimize function
foptim = @(p) logisticOptim(p,xx,yy);

options= optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true, 'HessianFcn', []);

params = fmincon(foptim, guess,[],[],[],[],lb,ub,[],options);  
%==========================================================================
%Xplot = linspace(min(xx), max(xx));
%plot(xx, yy, 'ob',Xplot, rlogistic4(params,Xplot),'-r',Xplot, rlogistic4(guess,Xplot),'-k');
%==========================================================================
%rescale params
params(1) = params(1 ) * max(spikes);
%==========================================================================
end

function [f,g,h] = logisticOptim(params, xx, yy)
% Calculate objective f

switch nargout
    case {0, 1}
        lf = nakarushton(params, xx);
    case 2
        [lf, lg] = nakarushton(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
    case 3
        [lf, lg, lH] = nakarushton(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
        h = (lg' * lg + reshape((lf-yy)' * lH(:,:), numel(params), numel(params)) )/size(xx, 1);
end

f = 0.5*(lf-yy)'*(lf-yy)/ size(xx, 1); %objective f



    
end