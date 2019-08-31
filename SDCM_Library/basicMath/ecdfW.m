%ABSTRACT
% Library function for SDCM. Weighted empirical CDF:
  function [SX,CDF,mass,SI] = ecdfW(X,W)
    %Initialize:
      bWeighted = nargin>=2 && all(size(W)==size(X));
      BValid = ~isnan(X);
      if(bWeighted) BValid = BValid & ~isnan(W); end
    %compute the CDF over sorted positions:
      [SX,SI] = sort(X(BValid));
      if(bWeighted)
        SW = W(BValid); SW = SW(SI);
        CDF = cumsum(abs(SW));
      else
        assert(nargin<2 || isscalar(W)&&W==1, 'W must be of size X or omitted (or ==1)');
        CDF = reshape(1:length(X), size(X));
      end
    %Normalize:
      if(~isempty(CDF))
        mass = CDF(end);
        CDF = CDF/mass;
      else
        mass = 0;
      end
  end
