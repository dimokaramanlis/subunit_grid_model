
% add repository path to MATLAB path
addpath(genpath('..\subunit_grid_model'))

%%
% set experiment folder for grating flashes
expfolder   = '20201112_252MEA_mouse_left_half_dorsal';

% load model data
modeltxt    = fullfile('..\subunit_grid_model', 'modelparameters', sprintf('%s.txt', expfolder));
cdata       = importdata(modeltxt, ' ', 1);
gridcenters = cdata.data(:,3:4);
mdlparams   = cdata.data(:, 5:end);
splitString = strsplit(cdata.textdata{1}, ' ');
textPart    = splitString{1};  % First part is the text
pxsize      = str2double(splitString{2});  % Second part is the number

% transform data to parameters for a single cell
icell       = 5;
mdlparams   = paramsToStructFlashes(mdlparams(icell, :), gridcenters(icell,:), 16/pxsize);
fprintf('Current cell is %d, id %d, model regularization used was %1.6f\n', ...
    icell, cdata.data(icell,1), cdata.data(icell,2))
