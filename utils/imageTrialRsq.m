function [ Rsqf, RsqTime, rfin, rTime] = imageTrialRsq( allTrialCounts )
%IMAGESELECTIVITYINDEX Calculates the selectivity to natural images
%   First defined in Quian Quiroga et al. (2007), also used by the Allen
%   Institue Brain Observatory.
%   Inputs:
%           allTrialCounts: n x Nt x p matrix, n is number of cells, Nt is number of trials
%   and p number of images
%   Outputs:
%           sIndex: n x 1 vector

[Ncells, Ntrials, Nframes] = size(allTrialCounts);
%==========================================================================
oddMean=reshape(mean(allTrialCounts(:,1:2:end,:),2), [Ncells Nframes]);
evenMean=reshape(mean(allTrialCounts(:,2:2:end,:),2), [Ncells Nframes]);

Rsq = 1-sum((oddMean-evenMean).^2,2)./sum((oddMean-repmat(mean(oddMean,2),[1 Nframes])).^2,2);
%Rsq = Rsq.*(Rsq>0);

Rsq2 = 1-sum((oddMean-evenMean).^2,2)./sum((evenMean-repmat(mean(evenMean,2),[1 Nframes])).^2,2);
%Rsq2 = Rsq2.*(Rsq2>0);

Rsqf = mean([Rsq Rsq2],2);

%==========================================================================
if nargout > 1
    earlyMean  = reshape(mean(allTrialCounts(:,1:floor(Ntrials/2),:),2), [Ncells Nframes]);
    lateMean   = reshape(mean(allTrialCounts(:,floor(Ntrials/2)+1:end,:),2), [Ncells Nframes]);

    Rsq  = 1-sum((earlyMean-lateMean).^2,2)./sum((earlyMean-repmat(mean(earlyMean,2),[1 Nframes])).^2,2);
    Rsq2 = 1-sum((earlyMean-lateMean).^2,2)./sum((lateMean-repmat(mean(lateMean,2),[1 Nframes])).^2,2);

    RsqTime = mean([Rsq Rsq2],2);
    
    if nargout > 2
        rTime = diag(corr(earlyMean', lateMean'));
        rfin  = diag(corr(oddMean', evenMean'));

    end

end
%==========================================================================
end

