%effect p values based on combined p values for correlations and its abs signal strength:
  function log10_p = pValue4SignalStrength(...
     effectAbsMean2D, nsBehindEffectAbsMean2D, effectAbsSD2D ...
    ,overallSignalAbsMean2D, nsBehindOverallSignalAbsMean2D, overallSignalAbsSD2D ...
  )
    %Probability of seeing this signal strength (via weighted 2-sample ttest):
      %Assume equal variance (we have the same global noise source, irrespective of the detected signatures's signal strength)
        dofs = nsBehindEffectAbsMean2D + nsBehindOverallSignalAbsMean2D - 2;
        combinedSD = sqrt((...
            (nsBehindEffectAbsMean2D-1) * effectAbsSD2D.^2 ...
          + (nsBehindOverallSignalAbsMean2D-1) * overallSignalAbsSD2D.^2 ...
        )/dofs);
        standardErrorOfTheMean = combinedSD .* sqrt(1./nsBehindEffectAbsMean2D + 1./nsBehindOverallSignalAbsMean2D);
      %t statistic:
        t = (effectAbsMean2D-overallSignalAbsMean2D) ./ standardErrorOfTheMean;
        %<-Note: the mean signal strength of the effect was calculated using the signs of correlations 
        %  to get a constructive sum, where the signal is consistent with the effect's correlations; 
        %  given perfect consistence, this equals the weightedMean(abs(signal)). Therefore we need to 
        %  use eDState.noiseEstimation.absMean2D and eDState.noiseEstimation.SD4absMean2D here for comparison. 
        %  (It does not make sense to try measuring the consistency by comparing to mean zero and 
        %  using eDState.noiseEstimation.noiseSD2D4currentL2Rs, since the algorithm always finds some correlations that 
        %  gives us a nonzero effect signal and this is biased. => Let the other p value for the correlations 
        %  and their underlying sample sizes handle significance testing of the consistency and compute here
        %  only the p value for the effect's signal strength conservatively, as if it would indeed 
        %  be weightedMean(abs(signal)).)
    %Compute p value:
      log10_p = log10(tcdf(double(t),double(dofs),'upper')); %upper tail as we are interested only in signals > abs noise niveau; others should get p->1.
      
  end
  
