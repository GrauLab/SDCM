%ABSTRACT
% Library function for SDCM. Estimate p values for weighted correlations via t statistics and
% optionally compute a single p value for multiple correlations via a KS test.

function [P4Rs, log10_p4Correlations] = pValues4Correlations(Rs, nsBehindRs, dim)
  %Estimate t statistics for the correlations:
    T4Rs = Rs .* sqrt(bsxfun(@times,max(0,nsBehindRs-2),1./(1-Rs.^2)));  
      %<-source with derivation from an uncorrelated bivariate normal distribution, see book "The analysis of physical measurements", Pugh Winslow, formula 12.93. 
      %<-other source: http://mathworld.wolfram.com/CorrelationCoefficientBivariateNormalDistribution.html
    %Special case nsBehindRs<=2 and Rs==1 (can happen in high NaN ratio scenarios): Here the above formula results in NaNs; define p value for correlations backed by only one sample as 1 here:
      T4Rs(nsBehindRs<=2) = 0;
    assert(isempty(T4Rs) || all(~isnan(T4Rs(:))), 'code validation: at least one T4Rs was NaN');

  %p values via tcdf:
    if(isscalar(nsBehindRs))
      P4Rs = 2*tcdf(-abs(T4Rs), max(eps,nsBehindRs-2)); %*2 to integrate both tails/signs (1 maps to 50:50 probability); max(eps,.) to prevent NaNs if nsBehindRs<=2 (can typically only happen in in high NaNs scenarios)
    else
      P4Rs = 2*bsxfun(@tcdf,-abs(T4Rs),max(eps,nsBehindRs-2)); %*2 to integrate both tails/signs (1 maps to 50:50 probability); max(eps,.) to prevent NaNs if nsBehindRs<=2 (can typically only happen in in high NaNs scenarios)
    end

  %combined p value via Kolmogorov Smirnov against T CDF:
    if(nargout>=2)
      if(dim==1) %nach unten aggregieren
        log10_p4Correlations = nan(1,size(Rs,2));
        for j = 1:size(Rs,2)
          %Due to different DOFs, Rs are not from the same t distribution. We need a linear combination of t statistics:
            Taxis = linspace(0,25,100);
            %approx by a single t distribution for the average DOFs:
              BInPerformanceSubspace = nsBehindRs(:,j)>1; %others were not correlated for performance resons.
              nu = mean(nsBehindRs(BInPerformanceSubspace,j));
              fcnNullCDF = griddedInterpolant(Taxis, 2*(tcdf(Taxis,nu)-0.5), 'linear', 'nearest'); %t distribution for t>=0; 2*(cdf-50%) utilizes symmetry (two-tailed).
              mint4comparison = tinv(1-0.5/2, nu); %for a more precise test, focus on the upper 50% of the null CDF (interesting genes are have high ts and are expected at the tail; any non-random deviations for lower ts are not interesting and should not lead to significance)
          log10_p4Correlations(j) = ksTestW(... %get one p for all Rs(:,j); does the detected signature have genes with significantly greater&more t values as expected by noise?
            ...test distribution and its weights:
                abs(T4Rs(BInPerformanceSubspace,j)), [] ...  %[] means unweighted, i.e. each t value has the same weight in the resulting CDF that is compared with fcnNullCDF.
            ...null distribution:
               ,fcnNullCDF, Inf ... %Inf means that the second input is considered a theoretically exact CDF, not just a sampling.
            ,mint4comparison ... %compare the upper tails of the distributions.
            ,'pSmallIfX1HasMoreHigherValues' ... %if T4Rs has more higher t-values than expected by the t distribution, correlations are significant.
          );
        end
      elseif(dim==2) %nach rechts aggregieren
        log10_p4Correlations = nan(size(Rs,1),1);
        for i = 1:size(Rs,1)
          %Due to different DOFs, Rs are not from the same t distribution. We need a linear combination of t statistics:
            Taxis = linspace(0,25,100);
            %approx by a single t distribution for the average DOFs:
              BInPerformanceSubspace = nsBehindRs(i,:)>1; %others were not correlated for performance resons.
              nu = mean(nsBehindRs(i,BInPerformanceSubspace));
              fcnNullCDF = griddedInterpolant(Taxis, 2*(tcdf(Taxis,nu)-0.5), 'linear', 'nearest'); %t distribution for t>=0; 2*(cdf-50%) utilizes symmetry (two-tailed).
              mint4comparison = tinv(1-0.5/2, nu); %for a more precise test, focus on the upper 50% of the null CDF (interesting genes are have high ts and are expected at the tail; any non-random deviations for lower ts are not interesting and should not lead to significance)
          log10_p4Correlations(i) = ksTestW(... %get one p for all Rs(i,:); does the detected signature have samples with significantly greater&more t values as expected by noise?
            ...test distribution and its weights:
                abs(T4Rs(i,BInPerformanceSubspace)), [] ...  %one p for all Rs of (points->their point on the effect curve)
            ...null distribution:
               ,fcnNullCDF, Inf ... %Inf means that the second input is considered a theoretically exact CDF, not just a sampling.
            ,mint4comparison ... %compare the upper tails of the distributions.
            ,'pSmallIfX1HasMoreHigherValues' ...
          );
        end
      end
    end
end
