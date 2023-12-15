function [ F, J ] = diffofgauss2dfun( params, xy)
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
theta = params(5); Ac    = params(6);
k     = params(7); As    = params(8);
aa    = (cos(theta)^2)/(2 * sx^2)  + (sin(theta)^2)/(2 * sy^2);
bb    = -sin(2*theta) /(4 * sx^2)  + (sin(2*theta))/(4 * sy^2);
cc    = (sin(theta)^2)/(2 * sx^2)  + (cos(theta)^2)/(2 * sy^2);
%==========================================================================

%==========================================================================
%Calculate function value
%==========================================================================
inExpC = -(aa * (xx - mx).^2 + 2 * bb * (xx - mx).*(yy - my) + cc * (yy - my).^2);
inExpS = inExpC * k^-2;
Fc     = Ac * exp(inExpC);
Fs     = As * exp(inExpS);
F      = Fc - Fs;
%==========================================================================

%==========================================================================
%Calculate Jacobian
%==========================================================================
if nargout>1
       

    dmx = 2 * (aa * (xx - mx) + bb * (yy - my)) .* (Fc - Fs * k^-2);
    dmy = 2 * (cc * (yy - my) + bb * (xx - mx)) .* (Fc - Fs * k^-2);
    
    dasx = - cos(theta)^2 * sx^-3;
    dbsx =   sin(2*theta) * sx^-3 / 2;
    dcsx = - sin(theta)^2 * sx^-3;
    
    dasy = - sin(theta)^2 * sy^-3;
    dbsy = - sin(2*theta) * sy^-3 / 2;
    dcsy = - cos(theta)^2 * sy^-3;
    
    dath = sin(2 * theta) * (-sx^-2 + sy^-2) / 2;
    dbth = cos(2 * theta) * (-sx^-2 + sy^-2) / 2;
    dcth = sin(2 * theta) * ( sx^-2 - sy^-2) / 2;
    
    dsx    = (Fc - Fs * (k^-2)).* ...
        (-dasx * (xx - mx).^2 - 2 * dbsx * (xx - mx).*  (yy - my) - dcsx * (yy - my).^2 );
    
    dsy    = (Fc - Fs * (k^-2)).* ...
        (-dasy * (xx - mx).^2 - 2 * dbsy * (xx - mx).*  (yy - my) - dcsy * (yy - my).^2 );
    dtheta = (Fc - Fs * (k^-2)) .* ...
        (-dath * (xx - mx).^2 - 2 * dbth * (xx - mx).*  (yy - my) - dcth * (yy - my).^2 );
    
    dk = 2 * Fs .* inExpC * k^-3;
    
    J = [dmx dmy dsx dsy dtheta exp(inExpC) dk -(exp(inExpS))];
    
    
end
%==========================================================================
end

