%ABSTRACT
% Library function for SDCM. Computes weighted means (faster version without NaN support).

  function [S,norm] = weightedMean_fastNoNans(X,W,dim,bIntegrate,bSupressWarnings)
    if(nargin<2) W=1; end
    if(nargin<3)
      if(size(X,2)==1) dim=1;
      elseif(size(X,1)==1) dim=2;
      else
        error('the summing dimension parameter is missing, but X is neither a row nor a column vector');
      end
    end
    if(nargin<4) bIntegrate=false; end
    if(nargin<5) bSupressWarnings=true; end %silently define 0/0=0 without a ton of warnings. 

    if(isscalar(W) && W==1)
      %Norm:
        if(nargout>=2) %only if required. 
          siz = size(X); siz(dim)=1;
          norm = size(X,dim)*ones(siz); 
        end
      if(bIntegrate)
        S = sum(X,dim);
      else
        S = mean(X,dim);
      end
    else
      assert(isscalar(W) || size(W,dim)==size(X,dim), 'weight vector is not aligned in aggregation dim');
      %Norm:
        if(~bIntegrate || nargout>=2) %only compute, if required. 
          norm = sum(abs(W),dim);
          if(any(norm(:)==0))
            if(~bSupressWarnings)
              warning('for %d/%d %s all weights were zero in the summing dimension %d; now adding epsilon to the norm, i.e. defining the mean of the empty set as zero to prevent NaNs'...
                ,sum(any(norm(:)==0)), numel(norm), iif(dim==1,'columns',iif(dim==2,'rows','vectors')), dim...
              );
            end
            norm = norm+eps;
          end
        end
      if(bIntegrate)
        S = sum(bsxfun(@times,W,X),dim);
      else
        S = bsxfun(@rdivide, sum(bsxfun(@times,W,X),dim), norm); %bsxfun instead of just ./ to also support nD input arrays.
      end
    end
    assert(~any(isnan(S(:))), 'code validation: NaNs detected at the end of weightedMean_fastNoNans');
  end
