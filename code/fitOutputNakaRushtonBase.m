function [fitprms,res] = fitOutputNakaRushtonBase(xvals, spikes, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
foptim = @(p) optimizeNonlinearity(p, xvals, spikes);
%==========================================================================
% setup guess from current params
if nargin < 3
    %guess = [max(spikes)/2 1 1];
    [svals, scents]  = getNonlinearity(xvals, spikes, 40, 1);
    guess            = fitNakaRushtonToSpikes(double(scents), double(svals));
else
    guess = varargin{1};
end
%==========================================================================
% paramerer bounds 
lb    = [  0    0             0   0];
ub    = [Inf  Inf max(xvals)*10 Inf];
%==========================================================================
% linear params
A = []; b = [];
%==========================================================================
options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients', false,...
    'HessianFcn', []);
[fitprms,res] = fmincon(foptim, guess, A, b, [], [], lb, ub, [], options);

%==========================================================================

end

function [f, g, H] = optimizeNonlinearity(p, xx, yy)


Np = sum(yy);

switch nargout
    case {0, 1}
        lf = nakarushtonbase(p, xx);% get spiking response
    case 2
        [lf, lg] = nakarushtonbase(p, xx);
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
    case 3
        [lf, lg, lH] = nakarushtonbase(p, xx);
        
        g =  -((yy./lf)' * lg - sum(lg, 1))/Np;
        
        ht1 = ((yy./(lf.^2)).*lg)' * lg;
        ht2 = reshape(- (yy./lf)' * lH(:, :), [numel(p), numel(p)]);
        ht3 = squeeze(sum(lH, 1));
        H = (ht1 + ht2 + ht3) / Np;
end

f = -(log(lf)'*yy - sum(lf))/Np; %objective f

end