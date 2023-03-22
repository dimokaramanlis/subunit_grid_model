function subwts = getGaussianSubWeights(params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
mx    = params.gaussparams(1);
my    = params.gaussparams(2);
sx    = params.gaussparams(3);
sy    = params.gaussparams(4);
theta = params.gaussparams(5);
aa = (cos(theta)^2)/(2 * sx^2)  + (sin(theta)^2)/(2 * sy^2);
bb = -sin(2*theta) /(4 * sx^2)  + (sin(2*theta))/(4 * sy^2);
cc = (sin(theta)^2)/(2 * sx^2)  + (cos(theta)^2)/(2 * sy^2);
%==========================================================================
xx = params.subcnts(:, 1);
yy = params.subcnts(:, 2);
%==========================================================================
subwts = exp(-(aa * (xx - mx).^2 + 2 * bb * (xx - mx).*(yy - my) + cc * (yy - my).^2));
%==========================================================================
end

