


function OutputMat = CelltoMatUE (inputcell, varargin)
%
% This function is an extension to MATLAB original cell2mat function but it
% works with cells with different dimentions. If the cell size is un-even
% it generate the zero in its place.
%
% ===============================Inputs====================================
%
%   inputcell : any cell.
%
%================================Output====================================
%
%   OutputMat : matrix of inputcell.
%
%
% written by Mohammad, 23.09.2013.
% added catch part on 5.06.2015
% added transpose part on 24.05.2016.

if size(inputcell,1) == 1
    inputcell = transpose(inputcell);
    transposeFlag = true;
else
    transposeFlag = false;
end;

maxLength=max(cellfun(@(x)numel(x),inputcell));
try
    OutputMat = cell2mat(cellfun(@(x)cat(2,x,nan(1,maxLength-length(x))),inputcell,'UniformOutput',false));
catch
    OutputMat = nan(length(inputcell),maxLength);
    for ii = 1: length(inputcell)
        OutputMat(ii,1:length(inputcell{ii})) = inputcell{ii};
    end;
end;

if transposeFlag
    OutputMat = transpose(OutputMat);
end;

end
