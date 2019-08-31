%ABSTRACT
% Library function for SDCM. Computes weighted means and supports NaNs.

function [S,norm] = weightedMean_supportsNaNs(X,W,dim,bIntegrate,bSupressWarnings)
  %Default parameters:
    if(nargin<2) W=1; end %unweighted, if no weights are provided.
    if(nargin<3)
      if(size(X,2)==1) dim=1;
      elseif(size(X,1)==1) dim=2;
      else
        error('the summing dimension parameter is missing, but X is neither a row nor a column vector');
      end
    end
    if(nargin<4) bIntegrate=false; end %if true, just integrate (do not devide by the norm).
    if(nargin<5) bSupressWarnings=true; end %silently define 0/0=0 without a ton of warnings. false; end

  if(isscalar(W) && W==1) %fast unweighted case:
    %Find NaNs:
      BNaNMask = isnan(X);
      X(BNaNMask) = 0;
    %Norm:
      if(~bIntegrate || nargout>=2) %only compute, if required. 
        norm = sum(~BNaNMask,dim); 
      end
    if(bIntegrate)
      S = sum(X,dim);
    else
      S = bsxfun(@rdivide, sum(X,dim), norm); %bsxfun instead of just ./ to also support nD input arrays.
    end
  else %weighted case:
    assert(isscalar(W) || size(W,dim)==size(X,dim), 'weight vector is not aligned in aggregation dim');
    %Find NaNs: support NaNs in X (are not summed by nansum, so their weights in the norm have to be zeroed, too):
      WX = bsxfun(@times, W, X);
      BNaN = isnan(WX);
      BAllNaNInSumDim = all(BNaN, dim);
      WX(BNaN) = 0;
    %Norm:
      if(~bIntegrate || nargout>=2) %only if required. 
        W(isnan(W)) = 0; %support NaNs in weights; needed since 0*NaN is still NaN.
        norm = sum(bsxfun(@times, double(~BNaN), abs(W)), dim); %support NaNs in X (are not summed by nansum, so their weights in the norm have to be zeroed, too).
        if(any(norm(:)==0 & ~BAllNaNInSumDim(:)))
          if(~bSupressWarnings)
            warning('for %d/%d %s all weights were zero in the summing dimension %d; now adding epsilon to the norm, i.e. defining the mean of the empty set as zero to prevent NaNs'...
              ,sum(any(norm(:)==0)), numel(norm), iif(dim==1,'columns',iif(dim==2,'rows','vectors')), dim...
            );
          end
          norm = norm+eps;
        end
      end
    if(bIntegrate)
      S = sum(WX, dim);
    else
      S = bsxfun(@rdivide, sum(WX, dim), norm); %bsxfun instead of just ./ to also support nD input arrays.
    end
    S(BAllNaNInSumDim) = NaN; %if all input in the summing dimension is NaN, so is their integral or mean.
  end
end
