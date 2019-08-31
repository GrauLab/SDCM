%ABSTRACT
% Library function for SDCM. Kolmogorov Smirnov test, supporting weighted distributions:

function [log10_p,ks] = ksTestW(X1,W1orMass1, X2orCDF2,W2orMass2, commonMinX,sTail, bPlot) %,sPSmallIfX1hasThisValuesRelativeToThisX2)
  %%.Initialization and parameter check:
    bUnweighted1 = isempty(W1orMass1);
    if(~bUnweighted1)
      W1 = W1orMass1;
      assert(all(size(W1)==size(X1)), 'W1 must be of size(X1) or [] for the equi-weighted case w_i==1');
      assert(all(W1(:)>=0 & W1(:)<=1), 'weights W1 must be in [0,1]');
    end
    bExternallyProvidedCDF2 = isa(X2orCDF2,'griddedInterpolant');
    if(~bExternallyProvidedCDF2)
      X2 = X2orCDF2;
      bUnweighted2 = isempty(W2orMass2);
      if(~bUnweighted2)
        W2 = W2orMass2;
        assert(all(size(W2)==size(X2)), 'W2 must be of size(X2)');
        assert(all(W2(:)>=0 & W2(:)<=1), 'weights W2 must be in [0,1]');
      end
    else
      fcnCDF2 = X2orCDF2;
      assert(isscalar(W2orMass2)&&W2orMass2>=0, 'if X2orCDF2 is the precomputed CDF, W2orMass2 must specify the mass behind it');
      mass2 = W2orMass2;
    end
    if(nargin<5 || isempty(commonMinX)) commonMinX = -Inf; end
    if(nargin<6 || isempty(sTail)) sTail = 'pSmallIfUnequal'; end
    if(nargin<7) bPlot = false; end

  %% Compute ks statistic:
    B1 = X1>=commonMinX; %also excludes NaNs.
      if(~any(B1(:)))
        %silently return p=1 for all-NaN inputs/warning('X1(X1>=commonMinX) was empty; returning p=1');
        log10_p=0; ks=0;
        return;
      end
    if(~bUnweighted1)
      [SX1,CDF1,mass1] = ecdfW(X1(B1),W1(B1));
    else
      [SX1,CDF1,mass1] = ecdfW(X1(B1));
    end
      if(~bExternallyProvidedCDF2) %performance/no need to interpolate if bExternallyProvidedCDF2, because we will use commonX=SX1 below.
        BDistinct1 = [true;diff(SX1(:))>eps(SX1(end))];
        fcnCDF1 = griddedInterpolant(SX1(BDistinct1), CDF1(BDistinct1), 'linear', 'nearest');
      end

    if(~bExternallyProvidedCDF2)
      B2 = X2>=commonMinX; %also excludes NaNs.
        if(~any(B2(:)))
          warning('X2(X2>=commonMinX) was empty; returning p=1');
          log10_p=0; ks=0;
          return;
        end
      if(~bUnweighted1)
        [SX2,CDF2,mass2] = ecdfW(X2(B2),W2(B2));
      else
        [SX2,CDF2,mass2] = ecdfW(X2(B2));
      end
        BDistinct2 = [true;diff(SX2(:))>eps(SX2(end))];
        fcnCDF2 = griddedInterpolant(SX2(BDistinct2), CDF2(BDistinct2), 'linear', 'nearest');
      commonX = unique([SX1(:);SX2(:)]);
    else
      commonX = SX1;
    end
    if(~bExternallyProvidedCDF2)
      cdf1 = fcnCDF1(commonX);
    else
      cdf1 = CDF1; %performance/commonX=SX1 in this case anyway.
    end
      %before the first point at min(SX1), define cdf1 as zero (instead of using nearest-neighbour interpolation):
        cdf1(commonX<=nanmin(SX1)) = 0; %use <= to get maximal ks if one CDF with few points start where the other is already high
        cdf1(commonX>=nanmax(SX1)) = 1;
    cdf2 = fcnCDF2(commonX);
      %before the first point at min(SX1), define cdf1 as zero (instead of using nearest-neighbour interpolation):
        if(~bExternallyProvidedCDF2)
          cdf2(commonX<=nanmin(SX2)) = 0; %use <= to get maximal ks if one CDF with few points start where the other is already high
          cdf2(commonX>=nanmax(SX2)) = 1;
        else
          cdf2(commonX<=min(fcnCDF2.GridVectors{1})) = 0; %use <= to get maximal ks if one CDF with few points start where the other is already high
          cdf2(commonX>=max(fcnCDF2.GridVectors{1})) = 1;
        end
      %if an externally provided CDF begins left of commonMinX, subtract the left mass and rescale, so that it starts with zero at commonMinX and ends with 100% (the comparison needs the same domain/support for both CDFs):
        if(bExternallyProvidedCDF2 && commonMinX>-Inf) 
          cdf2 = (cdf2-fcnCDF2(commonMinX))/(1-fcnCDF2(commonMinX));
        end
    switch(sTail)
      case 'pSmallIfX1HasMoreHigherValues'
        ks = max(-cdf1 + cdf2); %large if CDF1<CDF2 => more values(X1)>values(X2)
      case 'pSmallIfX2HasMoreHigherValues'
        ks = max(-cdf2 + cdf1); %large if CDF2<CDF1 => more values(X2)>values(X1)
      case 'pSmallIfUnequal'
        ks = abs(max(-cdf1 + cdf2));
    end

  %% Compute asymptotic p value for two sample-based CDFs (see kstest2.m)
    n1     =  mass1;
    n2     =  mass2;
    if(n2<inf)
      n      =  n1 * n2 /(n1 + n2);
    else
      n      =  n1; %if the second input is considered a theoretically exact CDF, n2->inf => n->n1.
    end
    lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * ks, 0);
    log10_p = -2*lambda.^2/log(10); % p = exp(-2 * lambda * lambda);
    if(isnan(log10_p))
      warning('code validation: isnan(log10_p) detected in ksTestW');
      keyboard;
    end

  %Simple plot if requested:
    if(bPlot)
      figure;
      plot(commonX,cdf1,'b.-'); 
      hold on; 
      plot(commonX,cdf2,'go-'); 
      axis tight;
      title(sprintf('ks=%0.3f, log_{10}(p)=%0.2f', ks,log10_p));
      %eDF({X1,X2},{W1,W2});
    end
end
