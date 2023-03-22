function gproj = boundedProjectedGradient(gorig, x, lb, ub)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

alpha = 1e-3;

p = x - alpha * gorig;

p(p < lb) = lb(p < lb);
p(p > ub) = ub(p > ub);

gproj = (x - p)/alpha;




end

