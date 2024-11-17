function nllpspk = neglogliperspike(spikes, preds)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nllpspk = -(log(preds(spikes>0))'*spikes(spikes>0) - sum(preds))/sum(spikes);

end

