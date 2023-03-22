function [ f, J, H] = rlogistic2( params, x )
%LOGISTICFUN evaluates modified logistic with params at x
%   Detailed explanation goes here
%==========================================================================
%define inputs
x     = x(:);
gamma = params(1); 
beta  = params(2);
%==========================================================================
%calculate function value
insignal = beta.*x + gamma;
eminus = exp(-insignal);
f = (1 + eminus).^-1;
%==========================================================================
if nargout>1
  
    J(numel(insignal), 2) = 0;
  
    J(:,1) = eminus .* (f.^2);
    J(:,2) = J(:,1).*x;
end
%==========================================================================
if nargout>2

    H = zeros(numel(x), numel(params),  numel(params));
    
    dnormpdf = J(:,1) .* (2 * eminus .* f - 1);
    dx = dnormpdf.*x;
    
    H(:,1,1) = dnormpdf;
    H(:,1,2) = dx;
    H(:,2,1) = dx;
    H(:,2,2) = dx .* x;
end


end

