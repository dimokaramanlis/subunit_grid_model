function [pts] = generateHexSubunitGrid(NsubMax)
%GENERATEHEXSUBUNITGRID Hexagonal grid of unit spacing centered at (0,0)
%   Detailed explanation goes here

RadMax = floor(max(roots([3 3 1-NsubMax])));

%Generate generic grid of NsubMax subunits
pts1 = squareGrid(RadMax * [-1 -1 1 1],           [0 0], [1 sqrt(3)]);
pts2 = squareGrid(RadMax * [-1 -1 1 1], [1/2 sqrt(3)/2], [1 sqrt(3)]);
pts = [pts1;pts2];

[~, isort] = sort(sqrt(sum(pts.^2,2)), 'ascend');
Nget = min(NsubMax, numel(isort));
pts = pts(isort(1:Nget), :);

end

