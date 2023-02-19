function [gauss] = getGaussFromNewParams(gparams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mx    = gparams(1); my    = gparams(2);
sx    = gparams(3); sy    = gparams(4);
theta = gparams(5);
aa = (cos(theta)^2)/(2 * sx^2)  + (sin(theta)^2)/(2 * sy^2);
bb = -sin(2*theta) /(4 * sx^2)  + (sin(2*theta))/(4 * sy^2);
cc = (sin(theta)^2)/(2 * sx^2)  + (cos(theta)^2)/(2 * sy^2);


S = inv(2*[aa bb; bb cc]);

gauss.mu    = [mx; my];
gauss.sigma = S ;


end

