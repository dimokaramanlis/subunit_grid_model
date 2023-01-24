function [ F, J ] = gauss2dfun( params, xy)
%GAUSSIAN2DFUN evaluates modified gaussian with params at x
%   Output: 
%       F: function value
%       J: Jacobian of the parameters
%Written by Dimos.

%==========================================================================
%Unwrap inputs
%==========================================================================
xx    = xy{1}(:);  yy = xy{2}(:);
mx    = params(1); my    = params(2); 
sx    = params(3); sy    = params(4);
theta = params(5); A     = params(6);
aa    = (cos(theta)^2)/(2 * sx^2)  + (sin(theta)^2)/(2 * sy^2);
bb    = -sin(2*theta) /(4 * sx^2)  + (sin(2*theta))/(4 * sy^2);
cc    = (sin(theta)^2)/(2 * sx^2)  + (cos(theta)^2)/(2 * sy^2);
%==========================================================================

%==========================================================================
%Calculate function value
%==========================================================================
inExp = -(aa * (xx - mx).^2 + 2 * bb * (xx - mx).*(yy - my) + cc * (yy - my).^2);
F     = A * exp(inExp);
%==========================================================================

%==========================================================================
%Calculate Jacobian
%==========================================================================
if nargout>1
       
    dmx = 2 * F .* (aa * (xx - mx) + bb * (yy - my));
    dmy = 2 * F .* (cc * (yy - my) + bb * (xx - mx));
    
    dasx = - cos(theta)^2 * sx^-3;
    dbsx =   sin(2*theta) * sx^-3 / 2;
    dcsx = - sin(theta)^2 * sx^-3;
    
    dasy = - sin(theta)^2 * sy^-3;
    dbsy = - sin(2*theta) * sy^-3 / 2;
    dcsy = - cos(theta)^2 * sy^-3;
    
    dath = sin(2 * theta) * (-sx^-2 + sy^-2) / 2;
    dbth = cos(2 * theta) * (-sx^-2 + sy^-2) / 2;
    dcth = sin(2 * theta) * ( sx^-2 - sy^-2) / 2;
    
    dsx    = F .* ...
        (-dasx * (xx - mx).^2 - 2 * dbsx * (xx - mx).*  (yy - my) - dcsx * (yy - my).^2 );
    dsy    = F .* ...
        (-dasy * (xx - mx).^2 - 2 * dbsy * (xx - mx).*  (yy - my) - dcsy * (yy - my).^2 );
    dtheta = F .* ...
        (-dath * (xx - mx).^2 - 2 * dbth * (xx - mx).*  (yy - my) - dcth * (yy - my).^2 );
    
    J = [dmx dmy dsx dsy dtheta exp(inExp)];
    
    
%     dmx=-(2*(1-rho^2))^-1 * (-2*(x-mx)/sx^2+2*rho.*(y-my)/(sx*sy));
%     J1=A*multi*exp(inExp).*dmx;
% 
%     dmy=-(2*(1-rho^2))^-1 * (-2*(y-my)/sy^2+2*rho.*(x-mx)/(sx*sy));
%     J2=A*multi*exp(inExp).*dmy;
% 
%     dsx1=-multi/sx;
%     dsx2=(-(2*(1-rho^2))^-1).*(-2*sx^-3*(x-mx).^2+2*sx^-2*rho*(x-mx).*(y-my)/sy);
%     J3=A*exp(inExp).*(dsx1+multi.*dsx2);
% 
%     dsy1=-multi/sy;
%     dsy2=(-(2*(1-rho^2))^-1).*(-2*sy^-3*(y-my).^2+2*sy^-2*rho*(x-mx).*(y-my)/sx);
%     J4=A*exp(inExp).*(dsy1+multi.*dsy2);
% 
%     dr1=multi*rho/(1-rho^2);
%     dr2=inExp*2*rho/(1-rho^2) + (-(2*(1-rho^2))^-1).*(- 2*(x-mx).*(y-my)/(sx*sy));
%     J5=A*exp(inExp).*(dr1+multi.*dr2);
% 
%     J6=multi*exp(inExp);
% 
%     J=[J1 J2 J3 J4 J5 J6];
end
%==========================================================================
end

