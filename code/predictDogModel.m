function [impreds, stadog] = predictDogModel(mdlparams, imEnsemble, nrx, nry)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ypix, xpix, Nims] = size(imEnsemble);

swt = mdlparams.surrwt / mdlparams.surrsc^2;

rfparams    = [mdlparams.gaussparams 1 mdlparams.surrsc swt];
% rfparams(1:2) = rfparams(1:2) + 0.5;
% centergauss = getGaussFromNewParams(rfparams);

% cel = getEllipse(centergauss, 4);
% %cel(2, :) = ypix - cel(2,:);
% maxLim    = floor(max(cel,[],2) + 0.5);
% minLim    = floor(min(cel,[],2) + 0.5);
% cutAround = max(maxLim - minLim);
% 
% rangeX = ceil((minLim(1)+maxLim(1)-cutAround)/2):ceil((minLim(1)+maxLim(1)+cutAround)/2);
% rangeX = rangeX(rangeX>0 & rangeX<=xpix);
% rangeY = ceil((minLim(2)+maxLim(2)-cutAround)/2):ceil((minLim(2)+maxLim(2)+cutAround)/2);
% rangeY = rangeY(rangeY>0 & rangeY<=ypix);
% 
% [cx, cy] = meshgrid(single(rangeX), single(rangeY));
[cx, cy] = meshgrid(single(nrx), single(nry));

stadog   = diffofgauss2dfun(rfparams, {cx(:), cy(:)});
stadog   = reshape(stadog,numel(nry), numel(nrx));
stadog   = stadog/sum(abs(stadog(:)));

impreds = stadog(:)' * reshape(imEnsemble(nry, nrx,:),...
    [numel(nry)* numel(nrx), Nims]);


end

