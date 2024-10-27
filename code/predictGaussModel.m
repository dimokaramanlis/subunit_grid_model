function [impreds, stagauss] = predictGaussModel(mdlparams, imEnsemble, nrx, nry)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[ypix, xpix, Nims] = size(imEnsemble);

imenuse = reshape(imEnsemble(nry, nrx, :), [numel(nry)*numel(nrx), Nims]);

rfparams    = [mdlparams.gaussparams 1];
%rfparams(1:2) = rfparams(1:2) + 0.5;

[cx, cy] = meshgrid(single(nrx), single(nry));
stagauss   = gauss2dfun(rfparams, {cx(:), cy(:)});
stagauss   = reshape(stagauss,numel(nry), numel(nrx));
stagauss   = stagauss/sum(abs(stagauss(:)));

impreds = stagauss(:)' * imenuse;



% [ypix, xpix, Nims] = size(imEnsemble);

% rfparams    = [mdlparams.gaussparams 1];
% rfparams(1:2) = rfparams(1:2) + 0.5;
% centergauss = getGaussFromNewParams(rfparams);
% 
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
% stadog   = gauss2dfun(rfparams, {cx(:), cy(:)});
% stadog   = reshape(stadog,numel(rangeY), numel(rangeX));
% stadog   = stadog/sum(abs(stadog(:)));
% 
% impreds = stadog(:)' * reshape(imEnsemble(rangeY, rangeX,:),...
%     [numel(rangeY)* numel(rangeX), Nims]);


end
