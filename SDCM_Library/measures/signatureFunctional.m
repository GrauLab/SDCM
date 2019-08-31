%ABSTRACT
% Library function for SDCM. Signature functional, i.e. a scalar score for 
% a signature. The score is high, if the interaction represented by the 
% signature is comprised of many and highly correlated genes and samples.


function score = signatureFunctional(signature)
  r = signature.signatureCorrInExtendedFocus./sqrt(1-signature.signatureCorrInExtendedFocus^2); %inpired by the t statistic corresponding to a correlation; better than just r as it results in a higher dynamic range of the score for rs near 1.
  s = sqrt(min(signature.signatureSizeByCorrSum4G,signature.signatureSizeByCorrSum4P)^0.9*max(signature.signatureSizeByCorrSum4G,signature.signatureSizeByCorrSum4P)^0.1); %sqrt of the signature "pixel area", but emphasis on the more narrow side (signatures involving many genes and many samples should get a higher score as signatures of the same area that are narrow in one order dimension, e.g. only affecting 4 samples...)
  score = r*s; %"correlation area integral"
  assert(~isnan(score), 'code validation: signature functional returned NaN');
end

