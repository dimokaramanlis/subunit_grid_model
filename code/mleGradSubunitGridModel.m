function [fout, insignal,gbatch] = mleGradSubunitGridModel(mdlparams, xvals, yy, ArgSubs, iuse,lambda,wd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargout <3
    Rf           = calcSubunitActivationAmplitudes(mdlparams, xvals(iuse, 1));
    subacts      = Rf .* ArgSubs(iuse, :);
    RsubsNlin    = rlogistic2gpu(mdlparams.subparams, subacts);
    RsubsNlin    = reshape(RsubsNlin, size(subacts));
    insignal     = RsubsNlin * mdlparams.subwts + mdlparams.bconst;
    %fall         =modifiedswish(mdlparams.outparams, insignal);
    fall          = nakarushton(mdlparams.outparams, insignal);
    fout = fall  + mdlparams.outbase;
end
%==========================================================================
if nargout > 2
    % forward and backward pass
    Ns = numel(mdlparams.subparams);
    outprms = mdlparams.outparams;
    [Rf, Rj]           = calcSubunitActivationAmplitudes(mdlparams, xvals(iuse, 1));
    subacts            = Rf .* ArgSubs(iuse, :);
    [RsubsNlin, jnlin] = rlogistic2gpu(mdlparams.subparams, subacts);
    RsubsNlin          = reshape(RsubsNlin, size(subacts));
    insignal           = RsubsNlin * mdlparams.subwts + mdlparams.bconst;
    %[fall, jf]         = modifiedswish(outprms, insignal);% get spiking response
    [fall, jf]          = nakarushton(outprms, insignal);
    fout = fall + mdlparams.outbase;

%     inexp = exp(-(outprms(3)*insignal + outprms(4)));
%     dfdz  = fall .* (outprms(2)./insignal + outprms(3)*inexp./(1+inexp));
    
    dfdz  = outprms(2) * fall .*(1-fall/outprms(1)) ./insignal;
    
    dNdg  = mdlparams.subparams(2) * RsubsNlin .* (1 - RsubsNlin);

    Jsigma = ((dNdg .*ArgSubs(iuse,:)) * mdlparams.subwts) .* Rj;
    
    Jsub  = permute(reshape(jnlin, [size(subacts) Ns]), [1 3 2]);
    Jsub  = reshape(Jsub, numel(iuse)*Ns, size(subacts,2)) * mdlparams.subwts;
    Jsub  = reshape(Jsub, numel(iuse), Ns);

    gall  = [dfdz.* [RsubsNlin Jsigma Jsub] jf ones(size(dfdz))];
    %==========================================================================
    Np = sum(yy);
    ycalc = yy(iuse);
    gbatch =  -((ycalc(ycalc>0)./fout(ycalc>0))' * gall(ycalc>0,:) - sum(gall, 1))/Np;
    Nw = numel(mdlparams.subwts);

    wuse  = wd*mdlparams.subwts;
    gbatch(1:Nw) = gbatch(1:Nw) +  2* lambda * wuse'.* ((mdlparams.subwts>0))';

    %gbatch(1:Nw) = gbatch(1:Nw) +  lambda .* ((mdlparams.subwts>0))';

    gbatch = gbatch';
    %==========================================================================
end
%==========================================================================
end