function [ f, J, H] = rlogistic3( params, x )
%LOGISTICFUN evaluates modified logistic with params at x
%   Detailed explanation goes here
%==========================================================================
%define inputs
x       = x(:);
gamma   = params(1); 
beta    = params(2); 
alpha   = params(3);
%==========================================================================
%calculate function value
fbef = (1 + exp(- beta * x  - gamma)).^-1;
f    = alpha * fbef;
%==========================================================================
if nargout>1
    J(numel(x), numel(params)) = 0;
    J(:,1) = f .* (1 - fbef);
    J(:,2) = J(:,1).*x;
    J(:,3) = f/alpha;
end
%==========================================================================
if nargout>2
    % get Hessian
    H(numel(x), numel(params),  numel(params)) = 0;        
    dnormpdf = f .* (2 * fbef.^2 + fbef + 1);
    
    dx = dnormpdf.*x;
    
    H(:,1,1) = dnormpdf;
    H(:,1,2) = dx;
    H(:,1,3) = J(:,1)/alpha;
    
    H(:,2,1) = dx;
    H(:,2,2) = dx .* x;
    H(:,2,3) = J(:,2)/alpha;

end
%==========================================================================
end

