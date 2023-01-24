function [params]= fitgaussrf(spx,spy,img, varargin)
%DOGRECEPTFIELD
%==========================================================================
%define ranges
cx = median(spx); 
cy = median(spy);
dx = abs(mean(diff(spx))); 
dy = abs(mean(diff(spy)));

[X,Y] = meshgrid((spx-cx)/dx,(spy-cy)/dy);
%==========================================================================
foptim=@(p) gaussOptim(p,{X(:), Y(:)},img(:));
%==========================================================================
%get guess
if nargin < 4
    guess = getGuess(X, Y, img);
else
    guessori = varargin{1};
    guess    = guessori;
    guess(1:2) = guess(1:2) - [cx, cy];
    guess(3:4) = guess(3:4)./ [dx dy];
    if numel(guess)< 6
        guess(6) = 1;
    end
end
%==========================================================================
if ~isfinite(foptim(guess))
    params = NaN(size(guess));
    return;
end
%==========================================================================
% if  ~isfinite(foptim(guess))
%     guess = [0 0 1 1 0 1];
% end
%==========================================================================
%Define lower and upper bounds for the parameters
rmax = max(range(X(:)), range(Y(:)))/2;
lb   = [guess(1) - rmax  guess(2) - rmax   0.1   0.1  -pi/4   0];
ub   = [guess(1) + rmax  guess(2) + rmax  rmax  rmax   pi/4 Inf];
%==========================================================================
%Setting up linear inequality constraints
A=[];
b=[]; 
%==========================================================================
%Perform the fit
options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
    'Display','off','SpecifyObjectiveGradient',true,'CheckGradients',false,...
    'FiniteDifferenceType','central');

[params, res] = fmincon(foptim, guess,A,b,[],[],lb,ub,[],options);
%==========================================================================
%Bring parameters to original scale
params(1)=params(1)*dx+cx; params(2)=params(2)*dy+cy;
params(3)=params(3)*dx; params(4)=params(4)*dy;
%==========================================================================
end


function [f,g] = gaussOptim(p,X,Y)
% Calculate objective f

[lf,lg] = gauss2dfun(p,X);

f = sum((lf-Y).^2)/2;

if nargout > 1 % gradient required
    g = (lf-Y)'*lg;
end
    
end

function guess = getGuess(X, Y, img)

[amp, im] = max(img(:));

rfinds=(img>3*std(img(:)));
if sum(rfinds)>0; guessmx=mean(X(rfinds)); guessmy=mean(Y(rfinds));
else, guessmx=X(im); guessmy=Y(im); end

guessArea=max([sum(rfinds(:)) 1]);
guessSigma=1.5*sqrt(guessArea)/2;

guessAc=(amp)/exp(1); 

guess=[guessmx guessmy guessSigma guessSigma 0 guessAc];

end