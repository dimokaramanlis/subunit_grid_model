function impreds = predictSingleSubunitModel(mdlparams, imEnsemble, nrx, nry)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[ypix, xpix, Nims] = size(imEnsemble);

imenuse = reshape(imEnsemble(nry, nrx, :), [numel(nry)*numel(nrx), Nims]);
[xuse, yuse]  = meshgrid(single(nrx)-0.5, single(nry)-0.5);

%iuse   = mdlparams.subwts>0;
iuse   = mdlparams.subwts>0;

subwts = mdlparams.subwts(iuse);
xcents = mdlparams.subcnts(iuse, 1);
ycents = mdlparams.subcnts(iuse,2);

% Get filters
pxarr = round(4*mdlparams.subsigma);
subfilters = (xuse(:) - xcents').^2 + (yuse(:) - ycents').^2;
subfilters(sqrt(subfilters) > pxarr) = NaN; % filter useless values
subfilters = exp(-subfilters/(2* mdlparams.subsigma^2));
if isfield(mdlparams, "subcentwt")
    cwt = mdlparams.subcentwt;
else
    cwt = 1;
end
swt        = mdlparams.subsurrwt/(mdlparams.subsurrsc^2)/cwt;
subfilters = subfilters - swt * subfilters.^(mdlparams.subsurrsc^-2);
subfilters(isnan(subfilters)) = 0;
%subfilters   = subfilters./sum(abs(subfilters),1);

subfilters  = subfilters./(2 * pi * mdlparams.subsigma^2);

% Do calculations
subacts  = cwt*subfilters' * imenuse;
%subacts  = subacts/max(subacts(:));

if numel(mdlparams.subparams) == 2
    subactsn = reshape(rlogistic2(mdlparams.subparams, subacts), size(subacts));
else
    subactsn = reshape(subnonlinbasis(mdlparams, subacts), size(subacts));
end

%impreds = subwts' * subactsn + mdlparams.bconst;

impreds = subwts' * subactsn;



% subfilters = (xuse(:) - xcents').^2 + (yuse(:) - ycents').^2;
% %subfilters(sqrt(subfilters) > pxarr) = NaN; % filter useless values
% subfiltersc = exp(-subfilters/(2* mdlparams.subsigma^2));
% subfilterss = subfiltersc.^(mdlparams.subsurrsc^-2);
% 
% 
% subacts  = (subfiltersc' * imenuse - ...
%     subfilterss' * imenuse*mdlparams.subsurrwt/(mdlparams.subsurrsc^2))/...
%     (2 * pi * mdlparams.subsigma^2);
% 
% subactsn = reshape(rlogistic2(mdlparams.subparams, subacts), size(subacts));
% impreds = subwts' * subactsn + mdlparams.bconst;
% 
% 
% 



        
end
