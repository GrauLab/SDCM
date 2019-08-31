%ABSTRACT
% Library function for SDCM. Compute the mask for dimensions included in subsequent computations
% based on the performance configuration (no point in including dimensions with noisy signal strengths
% or low correlations; speed up the computation of correlations in large datasets by removing them).

function BsPerformanceSubspace4Correlations = selectDimensionsInPerformanceSubspace(nG,nP, R4G,R4P, log10PNoise4signalStrengths4G,log10PNoise4signalStrengths4P, inInfo_searchStrategy)
  %Initialize:
    global eDState;
    sL = min(eDState.current.multiPassesLevel, inInfo_searchStrategy.nDefinedPasses); %config level.

  %Exclude all with too weak signal from the computation:
    threshold4G = log10(inInfo_searchStrategy.performance.alpha4Signal4inclusionInComputation4G(min(end,sL)));
    threshold4G = max(threshold4G, quantile(log10PNoise4signalStrengths4G, min(1,inInfo_searchStrategy.performance.minPerformanceSpaceSize4G(min(end,sL))/nG) ));
      BSufficientlyStrongSignal4G = log10PNoise4signalStrengths4G <= threshold4G;
    threshold4P = log10(inInfo_searchStrategy.performance.alpha4Signal4inclusionInComputation4P(min(end,sL)));
    threshold4P = max(threshold4P, quantile(log10PNoise4signalStrengths4P, min(1,inInfo_searchStrategy.performance.minPerformanceSpaceSize4P(min(end,sL))/nP) ));
      BSufficientlyStrongSignal4P = log10PNoise4signalStrengths4P <= threshold4P;
    
  %Rescue/Never exclude those with high absolute correlation to the current effect axes:
    BSufficientR4G = abs(R4G) >= inInfo_searchStrategy.performance.minAbsR4inclusionInCorrComputation(min(end,sL));
    BSufficientR4P = abs(R4P) >= inInfo_searchStrategy.performance.minAbsR4inclusionInCorrComputation(min(end,sL));
    
  %Combine: exclude noise, but never sufficiently (>0.5) correlated ones, even if they have a relatively weak signal:
    BsPerformanceSubspace4Correlations = {
      BSufficientlyStrongSignal4G | BSufficientR4G;
      BSufficientlyStrongSignal4P | BSufficientR4P;
    };
end

