function blockstimulus = returnFixMovie(screensize, imEnsemble, listfixations)
% returnFixMovie - Generates a movie of fixations using an image ensemble
%
% Syntax:  blockstimulus = returnFixMovie(screensize, imageEnsemble, listfixations)
%
% Inputs:
%    screensize - [Ny, Nx] size of the monitor, always [600, 800] for our experiments
%    imageEnsemble - ensemble of images, contains frozenImages and runningImages
%    listfixations - [3, Nframes] list of fixations (image index, x, y)
%                   for the frozen sequence it is frozenfixations, for the running sequence is runningfixations
%
% Outputs:
%    blockstimulus - [Ny, Nx, Nframes] movie of fixations, where the images presented at each frame are shifted based on gaze data
%
% Example:
%    blockstimulus = returnFixMovie([600, 800], imageEnsemble, listfixations)
%
% Other m-files required: none
% Subfunctions: getRanges
% MAT-files required: none
%--------------------------------------------------------------------------

Nx  = screensize(2);
Ny  = screensize(1);
[Nyim, Nxim, Nimages]  = size(imEnsemble);
%--------------------------------------------------------------------------
Nframes = size(listfixations, 2);
%--------------------------------------------------------------------------
blockstimulus = zeros([Ny,Nx, Nframes], 'single');

for ifix = 1:Nframes
    [xmin, xmax, ymin, ymax, xOr, yOr] = getRanges(...
        listfixations(2, ifix), listfixations(3, ifix),...
        Nx, Ny, Nxim, Nyim, 1:Nx, 1:Ny);
    blockstimulus(ymin:ymax, xmin:xmax, ifix) = imEnsemble(yOr, xOr, listfixations(1, ifix));
end
%--------------------------------------------------------------------------

end

function [xmin, xmax, ymin, ymax, xOr, yOr] = getRanges(trX, trY, Nxs, Nys, Nx, Ny, rx, ry)
     
ymin = (Nys / 2) - trY + 1; 
if ymin <= 1, ymin = 1; end

ymax = Ny + (Nys / 2) - trY; 
if (ymax > Nys),  ymax = Nys; end

xmin = (Nxs / 2) - trX + 1; 
if xmin <= 1, xmin = 1; end

xmax = Nx + (Nxs / 2) - trX; 
if (xmax > Nxs),  xmax = Nxs; end

% Rx    = xmin:xmax;
% rxuse = ismembc(rx, Rx);
% 
% Ry    = ymin:ymax;
% ryuse = ismembc(ry, Ry);

rxuse = (rx >= xmin & rx <=xmax);
ryuse = (ry >= ymin & ry <=ymax);

xOr = trX - (Nxs / 2) + rx(rxuse);
yOr = trY - (Nys / 2) + ry(ryuse);

xmin = find(rxuse, 1, 'first');
xmax = find(rxuse, 1, 'last');
ymin = find(ryuse, 1, 'first');
ymax = find(ryuse, 1, 'last');

end

