function [F, J] = funfitDoGModel(params, xvals, xorder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
gamma = params.outparams(1);
alpha = params.outparams(2);
%==========================================================================
% get gaussian prediction
[fc, fs, jc, js]  = calcDoGActivations(params, xvals);

kt_cent =  params.ktwts     * params.ktbasis';
kt_surr =  params.surrktwts * params.ktbasis';


allfc     = fc(xorder);
allfs     = fs(xorder);

allactiv = conv2(flip(kt_cent'), 1, allfc,'valid') + ...
    conv2(flip(kt_surr'), 1, allfs,'valid');
allactiv  = reshape(allactiv, [numel(allactiv), 1]);

cbaseproj = zeros(size(params.ktbasis,2), numel(allactiv));
sbaseproj = zeros(size(params.ktbasis,2), numel(allactiv));

for ii = 1:size(params.ktbasis,2)
    cbaseproj(ii,:) = reshape(conv2(flip(params.ktbasis(:,ii)), 1, allfc,'valid'),...
        [1,numel(allactiv)]);
    sbaseproj(ii,:) = reshape(conv2(flip(params.ktbasis(:,ii)), 1, allfs,'valid'),...
        [1,numel(allactiv)]);
end


% allfc     = fc(xorder);
% allfs     = fs(xorder);
% cbaseproj = (params.ktbasis' * allfc);
% sbaseproj = (params.ktbasis' * allfs);
% 
% allactiv    = params.ktwts * cbaseproj +  params.surrktwts * sbaseproj;
% allactiv    = allactiv';



alljc    = reshape(jc(xorder, :), [size(xorder,1), size(xorder,2)*6]);
alljs    = reshape(js(xorder, :), [size(xorder,1), size(xorder,2)*6]);

allj     = conv2(flip(kt_cent'), 1, alljc,'valid') + ...
    conv2(flip(kt_surr'), 1, alljs,'valid');
allj    = reshape(allj, [size(allj,1)*size(xorder,2), 6]);



% alljc    = reshape(jc(xorder, :), [size(xorder,1), size(xorder,2)*6]);
% alljs    = reshape(js(xorder, :), [size(xorder,1), size(xorder,2)*6]);
% allj    = kt_cent * alljc + kt_surr * alljs;
% allj    = reshape(allj, [size(xorder,2), 6]);

% get spiking response
fbef = (1 + exp(-(allactiv + gamma))).^-1;
F    = alpha * fbef;

%==========================================================================
if nargout > 1
    prefact = F.*(1-F/alpha);

    J1to6  = prefact .* allj;
    J      = [J1to6 prefact fbef prefact.* cbaseproj' prefact.* sbaseproj'];
end
%==========================================================================
end
