function mdlparams = modelToRAM(mdlparams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
allnames = fieldnames(mdlparams);
for ii = 1:numel(allnames)
    mdlparams.(allnames{ii}) = double(gather(mdlparams.(allnames{ii})));
end
end

