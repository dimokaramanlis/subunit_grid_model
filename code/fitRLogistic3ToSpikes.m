function params = fitRLogistic3ToSpikes(xvals, spikes)
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
    guess = [0 1 max(yy)];
else
    guess = [p(2) p(1) max(yy)];
end

lb = [-Inf -50   0];
ub = [Inf   50 Inf];
%==========================================================================
%optimize function
foptim = @(p) logisticOptim(p,xx,yy);

options= optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true, 'HessianFcn', [],...
    'CheckGradients', false);

params = fmincon(foptim, guess,[],[],[],[],lb,ub,[],options);  
%==========================================================================
%Xplot = linspace(min(xx), max(xx));
%plot(xx, yy, 'ob',Xplot, rlogistic3(params,Xplot),'-r',Xplot, rlogistic3(guess,Xplot),'-k');
%==========================================================================
%rescale params
params(3) = params(3) * max(spikes);
%==========================================================================
end

function [f,g,h] = logisticOptim(params, xx, yy)
% Calculate objective f

switch nargout
    case {0, 1}
        lf = rlogistic3(params, xx);
    case 2
        [lf, lg] = rlogistic3(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
    case 3
        [lf, lg, lH] = rlogistic3(params, xx);
        g = (lf-yy)' * lg /size(xx, 1);
        h = (lg' * lg + reshape((lf-yy)' * lH(:,:), numel(params), numel(params)) )/size(xx, 1);
end

f = 0.5*(lf-yy)'*(lf-yy)/ size(xx, 1); %objective f

end