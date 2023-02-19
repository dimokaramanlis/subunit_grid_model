function [ c ] = getEllipse( gauss, nsigma,detail)
%GETELLIPSE Summary of this function goes here
%   Detailed explanation goes here

if nargin<3; detail=5e1;end 
c = NaN(2,detail);
if sum(isnan(gauss.sigma))>0; return; end;
for j=1:size(c,2)
    alpha = 2*pi*(j-1)/(size(c,2)-1);
    c(:,j) = nsigma*sqrtm(gauss.sigma)*[cos(alpha);sin(alpha)]+gauss.mu; %do repmat to avoid for loop
end
c=[c c(:,1)];

% Converting input arguments into column vectors
%center=center(:)';
%sigmarule=sigmarule(:)';
%numpoints=ceil(numpoints);

% Calculates principal directions(PD) and variances (PV)
%[PD,PV]=eig(covmat);
%PV=diag(PV).^.5;

% Chooses points
%theta=linspace(0,2*pi,numpoints)';


% Construct ellipse
%elpt=[cos(theta),sin(theta)]*diag(PV)*PD';
%numsigma=length(sigmarule);
%elpt=repmat(elpt,1,numsigma).*repmat(sigmarule(floor(1:.5:numsigma+.5)),numpoints,1);
%elpt=elpt+repmat(center,numpoints,numsigma);

end

