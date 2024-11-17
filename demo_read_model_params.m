% This script loads and processes model parameters for a single cell from 
% a specified experiment. It reads the parameters from a text file in the repository, 
% extracts relevant information (grid centers, model parameters, pixel size), 
% and then converts the parameters into a structured format using either 
% 'paramsToStructFlashes' or 'paramsToStructFlicker' depending on the 
% type of stimulus used in the experiment ('gratingflashes' or 
% 'gratingflicker').

% add repository path to MATLAB path
addpath(genpath('..\subunit_grid_model'))

%%
% set experiment folder and select cell in experiment
expfolder   = '20220301_60MEA_marmoset_left_s1';
icell       = 5;

% load model data
modeltxt    = fullfile('..\subunit_grid_model', 'modelparameters', sprintf('%s.txt', expfolder));
cdata       = importdata(modeltxt, ' ', 1);
gridcenters = cdata.data(:,3:4);
mdlparams   = cdata.data(:, 5:end);
splitString = strsplit(cdata.textdata{1}, ' ');
textPart    = splitString{1};              % First part is the text
pxsize      = str2double(splitString{2});  % Second part is the number

% transform data to parameters for a single cell
switch textPart
    case 'gratingflashes'
        mdlparams   = paramsToStructFlashes(mdlparams(icell, :), gridcenters(icell,:), 16/pxsize);
    case 'gratingflicker'
        gflickerdata = load(fullfile('..\subunit_grid_model', expfolder, 'gratingflicker_data.mat'));
        mdlparams    = paramsToStructFlicker(...
            mdlparams(icell, :), gridcenters(icell,:), gflickerdata.ktbas, 16/pxsize, 7.5/pxsize);
end

fprintf('Current cell is %d, id %d, model regularization used was %1.6f\n', ...
    icell, cdata.data(icell,1), cdata.data(icell,2))

