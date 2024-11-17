function [pts] = generateHexSubunitGrid(NsubMax)
%GENERATEHEXSUBUNITGRID Generates a hexagonal grid of subunit locations.
%
%   This function creates a hexagonal grid of points with unit spacing, 
%   centered at (0,0). The grid is designed to accommodate up to NsubMax 
%   subunits, ensuring that the points are arranged in a compact hexagonal 
%   pattern. 
%
%   Inputs:
%       NsubMax - The maximum number of subunits to include in the grid.
%
%   Outputs:
%       pts     - An N x 2 matrix of (x, y) coordinates for the subunit 
%                 locations in the hexagonal grid. N is the number of 
%                 subunits, which may be less than or equal to NsubMax.
%
%   Explanation:
%       The function first determines the necessary radius of the hexagonal 
%       grid to accommodate NsubMax subunits. It then generates two 
%       square grids with appropriate offsets to create the hexagonal 
%       pattern. The points are sorted by their distance from the origin 
%       to ensure a compact arrangement. Finally, the function returns 
%       the coordinates of the NsubMax closest points.

RadMax = floor(max(roots([3 3 1-NsubMax])));

%Generate generic grid of NsubMax subunits
pts1 = squareGrid(RadMax * [-1 -1 1 1],           [0 0], [1 sqrt(3)]);
pts2 = squareGrid(RadMax * [-1 -1 1 1], [1/2 sqrt(3)/2], [1 sqrt(3)]);
pts = [pts1;pts2];

[~, isort] = sort(sqrt(sum(pts.^2,2)), 'ascend');
Nget = min(NsubMax, numel(isort));
pts = pts(isort(1:Nget), :);

end

