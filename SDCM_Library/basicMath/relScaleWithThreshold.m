%ABSTRACT
% Library function for SDCM. Transforms X into relative values and supports a threshold.
  function relX = relScaleWithThreshold(X, relXthreshold, relX4one, dim)
    if(nargin<4)
      if(isrow(X)) dim=2; elseif(iscolumn(X)) dim=1; else error('if relScaleWithThreshold is not called with a 1D vector, the dim argument is required'); end
    end
    relX = bsxfun(@times,X,1./(max(abs(X),[],dim)+eps)); %+eps to prevent 0/0=NaN for exact zero vectors.
    relX(bsxfun(@lt, abs(relX), relXthreshold)) = 0; %bsxfun to support multiple relXthreshold along the dual dim.
    relX = sign(relX) .* min(1, bsxfun(@rdivide,abs(relX),relX4one));
  end
