function params = fitRLogisticToSpikes(xvals, spikes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
%unwrap and standardize variables
spikes = spikes(:); 
xx     = xvals(:);
yy     = spikes/max(spikes);
%==========================================================================
yylog  = log((yy+1e-2)./(max(yy)-yy+1e-2));
p      = polyfit(xx, yylog, 1); %linear fit to get the estimates for the guess

if any(isnan(p))
    guess = [0 0 1 max(yy)];
else
    guess = [0 p(2) p(1) max(yy)];
end

lb = [-Inf -Inf -Inf   0];
ub = [Inf  Inf  Inf Inf];
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
params([1 4]) = params([1 4]) * max(spikes);
%==========================================================================
end

function [f,g,h] = logisticOptim(params, xx, yy)
% Calculate objective f

switch nargout
    case {0, 1}
        lf = rlogistic4(params, xx);
    case 2
        [lf, lg] = rlogistic4(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
    case 3
        [lf, lg, lH] = rlogistic4(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
        h = (lg' * lg + reshape((lf-yy)' * lH(:,:), numel(params), numel(params)) )/size(xx, 1);
end

f = 0.5*(lf-yy)'*(lf-yy)/ size(xx, 1); %objective f



    
end