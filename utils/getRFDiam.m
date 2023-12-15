function [ diamum ] = getRFDiam( gauss, nsigma, pixelsize)
%GETRFDIAM Summary of this function goes here
%   Input:
%       gauss: should be in screen units

if isempty(gauss); diamum=NaN;  
else; diamum=2*nsigma*pixelsize*det(gauss.sigma)^(1/4); end

end
 
