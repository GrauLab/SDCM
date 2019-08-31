%ABSTRACT
% Library function for SDCM. Return a structure of basic math 
% functions. Use fast funcitons without NaN support, if no
% NaN support is required. Otherwise return nan-robust versions.

function fcns = getBasicMathFunctions(bSupportNaNs)
  fcns = struct();
  if(bSupportNaNs)
    fcns.sum = @nansum;
    fcns.mean = @nanmean;
    fcns.meanW = @weightedMean_supportsNaNs;
    fcns.var = @(X,dim)nanvar(X,0,dim);
  else %use nonNan verions if possible for performance reasons:
    fcns.sum = @sum;
    fcns.mean = @mean;
    fcns.meanW = @weightedMean_fastNoNans;        
    fcns.var = @(X,dim)var(X,0,dim);
  end
  fcns.uncenteredVar = @(X,dim)fcns.mean((X-0).^2,dim); %Note: this is equivalent to the POWER of the signal.
  fcns.uncenteredVarW = @(X,W,dim)fcns.meanW((X-0).^2,W,dim); %Note: this is equivalent to the weighted POWER of the signal.
  fcns.euclidW = @(X,W,dim)fcns.meanW(X.^2,W.^2,dim,true).^(1/2);
end
