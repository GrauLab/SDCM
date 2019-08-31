%ABSTRACT
% Library function for SDCM. Signature focus, i.e. weights that are specific 
% and representative for the top genes/samples in a detected signature.
% Initial weights are based on the relative signal strengths. As soon as 
% correlations are availble, r*(1-p)� are used as more specific base weights.
% For the signature focus, base weights >=50% are already mapped to full 
% membership (w=1). For signature size estimation and regression weights, 
% this is not optimal; herefore, the extended signature focus keeps the full
% dynamic range (only base weights of one yield w=1) in the default config.

function Ws = calcSignatureFocus(signal, Ps4signal, Rs, Ps4Rs, dim, focusingConfig, sL)
  %initial weights based on signal strengths (step 1, while no correlations are available yet):
    if(sum(size(Rs))==0)
      %assert(~isempty(signal), 'signatureFocus.step1 needs signal to compute Ws'); 
      Ws = signedPower(relScaleWithThreshold(...
         signal ... %<-Note: we do not use Ps4signal as the true noise level is unknown; hence, Ps4signal (especially for k==1) are not reliable for all types of signals (keep the focus formula simple / do not rely on noise estimation for initial weights). We only need to focus roughly here anyway (these are only initial weights; the refined correlation-based weights are computed as defined below for every candidate, as soon as correlations have been computed in precomputeCandidatesInParallel).
        ,focusingConfig.relSignal.threshold(min(end,sL))...
        ,focusingConfig.relSignal.sufficient4fullWeight(min(end,sL))...
      ,dim),focusingConfig.relSignal.exponent(min(end,sL)));
      %map NaN weights as zero weights (for NaNs in the signal):
        Ws(isnan(Ws)) = 0;
    end

  %standard weights based on correlations and their p values:
    if(sum(size(Rs))>0)
      %Base weights:
        if(~focusingConfig.relCorr.bEnabled) %dev.note/reproducability feedback: the refined focus uses correlations relative to max(abs(Rs)), but the extended focus deoes not rescale correlations.
          Ws = Rs.*(1-Ps4Rs).^2; %base weights are correlations times a factor that is near one for significant correlations.
        else
          Ws = signedPower(relScaleWithThreshold(...
             Rs.*(1-Ps4Rs).^2 ... %base weights are correlations times a factor that is near one for significant correlations.
            ,focusingConfig.relCorr.threshold(min(end,sL))... %default=0
            ,focusingConfig.relCorr.sufficient4fullWeight(min(end,sL))... %default=0.5 for the signature focus, 1 for the extended signature focus
          ,dim),focusingConfig.relCorr.exponent(min(end,sL))); %default=1
        end
      %exclude noise dimensions (prevent a dominance of the noisy masses in high-dim datasets), using a gradual cutoff at 67% of the quantile axis:
        [SabsWs,SI] = sort(abs(Ws),dim);
        quantileAxis = linspace(1/size(Ws,dim),1,size(Ws,dim))';
          if(dim==2) quantileAxis = quantileAxis'; end
        BValid = bsxfun(@gt,SabsWs,focusingConfig.relCorr.cutNoiseBelowQuantileAxisRatio(min(end,sL))*quantileAxis); %default .cutNoiseBelowQuantileAxisRatio=0.67 for the signature focus, 0.4 for the extended signature focus.
          if(dim==1)
            BValid(sub2ind(size(Ws), SI, repmat(1:size(Ws,2),size(Ws,1),1))) = BValid;
          else
            BValid(sub2ind(size(Ws), repmat((1:size(Ws,1))',1,size(Ws,2)), SI)) = BValid;
          end
        Ws = Ws.*BValid;
      %Checks:
        assert(all(size(Ws)==size(Rs)), 'code validation: size(Ws)~=size(Rs)');
        assert(all(~isnan(Ws(:))), 'code validation: some weights were NaN'); %correlation-based weights cannot be NaN.
    end
    
  %Handle special cases:
    %we need at least 2 points for interpolation/bi-monotonic regression; otherwise zero the focus:
      BValid = sum(Ws~=0, dim)>=2;
      if(dim==1)
        Ws(:,~BValid) = 0;
      else
        Ws(~BValid,:) = 0;
      end
    %Warn if we have a zero focus everywhere:
      if(  focusingConfig.bWarnOnZeroWeights(min(end,sL)) ...
        ...&& any(sum(abs(Ws)==0,dim) >= size(Ws,dim)-1) ...
        && sum(sum(abs(Ws)==0,dim) >= size(Ws,dim)-1)/prod(illa(size(Ws),dim,1))>=0.05 ... %only warn, if this happens for >=5%
      )
        warning('Weights for %d/%d %s were all zero (or all but one), i.e. the current focusing or performance setting essentially excludes them. If this warning occurs during STEP 2 (signature focusing) and there are still sufficient member candidates to let the generalization converge, no information is lost and this warning can be safely ignored. Otherwise consider less tight settings or disable the performance subspace entirely; also check that your input data does not contain all zero or all nan rows or columns.'...
          ,sum(sum(abs(Ws)==0,dim) >= size(Ws,dim)-1)...
          ,prod(illa(size(Ws),dim,1)) ...
          ,iif(dim==1,'columns','rows')...
        );
      end
end

