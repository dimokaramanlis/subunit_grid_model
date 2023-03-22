function outplot = linedogplot(sigmac, surrsc, surrwt, xvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

swt = surrwt / surrsc^2;
signalc = exp(-xvals.^2/(2 * sigmac^2));
signals = exp(-xvals.^2/(2 * (surrsc * sigmac)^2));
outplot = signalc - swt * signals;


end

