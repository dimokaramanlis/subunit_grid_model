function c = getEllipseFromNewParams(gparams, nsigma, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3; detail=5e1;
else, detail = varargin{1}; end 

mx    = gparams(1); my    = gparams(2);
sx    = gparams(3); sy    = gparams(4);
theta = gparams(5);
aa = (cos(theta)^2)/(2 * sx^2)  + (sin(theta)^2)/(2 * sy^2);
bb = -sin(2*theta) /(4 * sx^2)  + (sin(2*theta))/(4 * sy^2);
cc = (sin(theta)^2)/(2 * sx^2)  + (cos(theta)^2)/(2 * sy^2);

if any(isnan([aa bb; bb cc]), 'all')
    S = NaN(2, 2);
else
    S = inv(2*[aa bb; bb cc]);
end
    
c = NaN(2,detail-1);
if sum(isnan(S))>0; return; end
for j=1:size(c,2)
    alpha = 2*pi*(j-1)/(size(c,2)-1);
    c(:,j) = nsigma*sqrtm(S)*[cos(alpha);sin(alpha)]+[mx;my]; %do repmat to avoid for loop
end
c=[c c(:,1)];
end

