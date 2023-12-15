function fitprms = fitOutputNonlinearityML4(xvals, spikes, guess)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
foptim = @(p) optimizeNonlinearity(p, xvals, spikes);
%==========================================================================
% setup guess from current params

%==========================================================================
% paramerer bounds 
lb    = [-Inf -Inf    0    0];
ub    = [Inf   Inf  Inf  Inf];
%==========================================================================
% linear params
A = []; b = [];
%==========================================================================
if guess(1)< 0
    guess(1)=0;
end
%==========================================================================
options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients', false,...
    'HessianFcn', 'objective', 'MaxIter', 200);
fitprms = fmincon(foptim, guess, A, b, [], [], lb, ub, [], options);

%==========================================================================

end

function [f, g, H] = optimizeNonlinearity(p, xx, yy)


Np = sum(yy);

switch nargout
    case {0, 1}
        lf = rlogistic4(p, xx);% get spiking response
    case 2
        [lf, lg] = rlogistic4(p, xx);
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
    case 3
        [lf, lg, lH] = rlogistic4(p, xx);
        
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
        
        ht1 = ((yy./(lf.^2)).*lg)' * lg;
        ht2 = reshape(- (yy./lf)' * lH(:, :), [numel(p), numel(p)]);
        ht3 = squeeze(sum(lH, 1));
        H = (ht1 + ht2 + ht3) / Np;
end

f = -(log(lf)'*yy - sum(lf))/Np; %objective f

end