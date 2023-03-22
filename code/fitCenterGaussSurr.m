function [gparams,rftofit] = fitCenterGaussSurr(mdlparams)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


pxranges = round(mdlparams.gridcent) + (-mdlparams.cwindow:mdlparams.cwindow)';
[xuse, yuse] = meshgrid(pxranges(:,1), pxranges(:,2));


subfilters = (xuse(:)-mdlparams.subcnts(:,1)').^2 +...
    (yuse(:)-mdlparams.subcnts(:,2)').^2;
subfilters = exp(-subfilters/(2*mdlparams.subsigma^2));

rftofit = reshape(subfilters*mdlparams.subwts, size(xuse));
rftofit = imgaussfilt(rftofit, 4*7.5/mdlparams.pxsize);
rftofit = rftofit - min(rftofit(:));
rftofit = rftofit/max(rftofit(:));
gparams = fitgaussrf(pxranges(:,1), pxranges(:,2), double(gather(rftofit)));


end

