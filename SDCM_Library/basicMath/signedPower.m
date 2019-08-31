%ABSTRACT
% Library function for SDCM. Computes sign(X).*|X|.^exponent.

  function X = signedPower(X,exponent)
    if(exponent==1)
      return; %performance.
    else
      X = sign(X).*abs(X).^exponent;
    end
  end
