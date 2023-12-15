function [F, J] = funfitGaussModel(params, xvals, xorder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
gamma = params.outparams(1);
alpha = params.outparams(2);
%==========================================================================
% get gaussian prediction
[fin, jin]  = calcGaussianActivationsGrating(params.gaussparams, xvals);

kt = params.ktwts * params.ktbasis';

allfin      = fin(xorder);
allactiv = conv2(flip(kt'), 1, allfin,'valid');
allactiv  = reshape(allactiv, [numel(allactiv), 1]);

baseproj = zeros(size(params.ktbasis,2), numel(allactiv));
for ii = 1:size(params.ktbasis,2)
    baseproj(ii,:) = reshape(conv2(flip(params.ktbasis(:,ii)), 1, allfin,'valid'),...
        [1,numel(allactiv)]);
end

% tic;
% Nt = numel(kt);
% baseproj2 = conv2(reshape(allfin,[numel(allfin),1]),flip(params.ktbasis,1),'full');
% baseproj2 = reshape(baseproj2(1:end-Nt+1,:), [size(xorder), size(params.ktbasis,2)]);
% baseproj2 = reshape(baseproj2(Nt:end, :, :), [numel(allactiv), size(params.ktbasis,2)])';
% toc;

% allfin      = fin(xorder);
% baseproj    = (params.ktbasis' * allfin);
% allactiv    = params.ktwts * baseproj;
% allactiv    = allactiv';
% 
alljs    = reshape(jin(xorder, :), [size(xorder,1), size(xorder,2)*5]);
alljs    = conv2(flip(kt'), 1, alljs,'valid');
alljs    = reshape(alljs, [size(alljs,1)*size(xorder,2), 5]);

% alljs    = reshape(jin(xorder, :), [size(xorder,1), size(xorder,2)*5]);
% alljs    = reshape(kt * alljs, [size(xorder,2), 5]);

% get spiking response
fbef = (1 + exp(-(allactiv + gamma))).^-1;
F    = alpha * fbef;

%==========================================================================
if nargout > 1
    prefact = F.*(1-F/alpha);

    J1to5  = prefact .* alljs;
    J      = [J1to5 prefact fbef prefact.* baseproj'];
end
%==========================================================================
end
