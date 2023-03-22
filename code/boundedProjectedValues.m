function x = boundedProjectedValues(x, lb, ub)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


x(x < lb) = lb(x < lb);
x(x > ub) = ub(x > ub);

end

