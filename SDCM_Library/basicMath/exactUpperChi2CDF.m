function ln_p = exactUpperChi2CDF(chi2,nu)
  %Special cases and quick return for scalar inputs:
    bIsScalar = isscalar(chi2);
    BIsInf = isinf(chi2);
    if(bIsScalar && BIsInf && nu>0) %if it already underflowed on caller level
      ln_p = -Inf;
      return;
    end
    BIsZero = chi2==0;
    if(bIsScalar && BIsZero && nu>0)  %if p==1 based on the combination of an empty set of p values (sum(Ws)==0 on caller level)
      ln_p = 0;
      return;
    end
  %Use scaled high precision integrator:
    %    p = chi2cdf(chi,nu,'upper');
    %<=> p = gamcdf(chi,nu/2,2,'upper'); %cf. chi2cdf.m
    %<=> p = gammainc(x, a,'upper') with x:=chi4g/(nu4g/2), a:=nu4g/2, b:=2 %cf. gamcdf.m
    %<=> p = gammainc(x, a,'scaledupper')/(gamma(a+1)*exp(x)/x^a) with x:=chi4g/(nu4g/2), a:=2 %cf. doc gammainc.
    %<- note: gammainc(x, a,'scaledupper') does not include the factor gamma(a+1)*exp(x)/x^a to prevent underflow to zero.
    a = nu/2; b = 2;
    x = chi2/b;
    ln_p = log(gammainc(x, a,'upper'));
    bUseScaledVersion = ln_p < -100;
    if(bUseScaledVersion)
      ln_scaled_p = log(gammainc(x, a,'scaledupper')); %includes factor Gamma(a+1)*exp(x)/x^a to prevent underflow to zero.
      ln_scaleFactor = gammaln(a+1) + x - a.*log(x);
      ln_p = ln_scaled_p - ln_scaleFactor;
    end
  %Validation:
    BIsNaN = isnan(ln_p);
    if(any(BIsNaN(:)))
      warning('code validation: any(isnan(ln_p(:))) detected, double underflow or allNaNs input?');
    end
end
