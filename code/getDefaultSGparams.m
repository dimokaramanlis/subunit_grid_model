function opts = getDefaultSGparams(modeltype)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


opts = struct();

opts.beta1     = .9;
opts.beta2     = .999;
opts.epsmall   = 1e-6;
opts.showfig   = false; % whether to show fit progress
switch modeltype
    case 'flicker'
        opts.lambda    = 50e-5;
        opts.batchsize = 2000;
        opts.eta       = 0.02;
        opts.Nepochs   = 40;
    case 'flashes'
end

end

