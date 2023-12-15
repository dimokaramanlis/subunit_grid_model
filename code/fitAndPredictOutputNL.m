function preds = fitAndPredictOutputNL(gensfit, genspred, spikesfit, funtype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[svals, scents]  = getNonlinearity(gensfit, spikesfit, 40, 1);

switch funtype
    case 'logistic4'
        guess       = fitRLogisticToSpikes(double(scents), double(svals));
        nlnparams   = fitOutputNonlinearityML4(double(gensfit), double(spikesfit),  ...
            [0 guess(2:end)]);
        funoutput   = @(x) rlogistic4(nlnparams, x);
    case 'nakarushton'
        
end

preds         = funoutput(genspred);



end

