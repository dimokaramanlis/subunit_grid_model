function [h_plot] = plotLogTuning2D(xvals, yvals, tuningcurve)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


yy       = zeros(numel(yvals)+1, 1);
alldiffsy = diff(yvals);
yy(2:end-1, 1) = yvals(1:end-1) + alldiffsy/2;
yy(1)   = yvals(1) - alldiffsy(1)/2;
yy(end) = yvals(end) + alldiffsy(end)/2;

xx       = zeros(numel(xvals)+1, 1);
alldiffsx = diff(xvals);
xx(2:end-1, 1) = xvals(1:end-1) + alldiffsx/2;
xx(1)   = xvals(1) - alldiffsx(1)/2;
xx(end) = xvals(end) + alldiffsx(end)/2;


xall = zeros(4, numel(tuningcurve));
yall = zeros(4, numel(tuningcurve));

idx = 1;
for iy = 1:numel(yvals)
    ycoords = yy(iy:iy+1);
    for ix = 1:numel(xvals)
        xcoords = xx(ix:ix+1);
        [xadd, yadd] = meshgrid(xcoords, ycoords);
        xall(:, idx) = xadd(:);
        yadd(:,2) = flip(yadd(:,2));
        yall(:, idx) = yadd(:);
        idx = idx + 1;
    end
end


vrs = [xall(:) yall(:)];
fcs = reshape(1:size(vrs,1), [4, numel(tuningcurve)])';
cdatause = tuningcurve';

ax = gca; ax.XScale = 'log'; 
ylim([min(yy) max(yy)]); xlim([min(xx) max(xx)]);
patch('Faces', fcs, 'Vertices', vrs, 'FaceColor', 'flat', ...
    'CData', cdatause(:), 'EdgeColor', 'none');


end

