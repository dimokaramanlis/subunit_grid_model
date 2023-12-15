function [ f, J, H] = rlogistic4( params, x )
%LOGISTICFUN evaluates modified logistic with params at x
%   Detailed explanation goes here
%==========================================================================
%define inputs
x       = x(:);
epsilon = params(1); 
gamma   = params(2); 
beta    = params(3); 
alpha   = params(4);
%==========================================================================
%calculate function value
insignal = beta.*x + gamma;
eminus   = exp(-insignal);
f = epsilon + alpha./(1 + eminus);
%==========================================================================
if nargout>1
    J = ones(numel(x), numel(params));
    %J1 = ones(size(x));
    J(:,2) = eminus .* (f - epsilon).^2 / alpha;
    J(:,3) = J(:,2).*x;
    J(:,4) = (f - epsilon)/ alpha;
end
%==========================================================================
if nargout>2
    % get Hessian
    H = zeros(numel(x), numel(params),  numel(params));    
    dnormpdf = J(:,2) .* (2 * eminus .* J(:,4) - 1);
    dx = dnormpdf.*x;
    
    H(:,2,2) = dnormpdf;
    H(:,2,3) = dx;
    H(:,2,4) = J(:,2)/alpha;
    
    H(:,3,2) = dx;
    H(:,3,3) = dx .* x;
    H(:,3,4) = J(:,3)/alpha;

end
%==========================================================================
end

