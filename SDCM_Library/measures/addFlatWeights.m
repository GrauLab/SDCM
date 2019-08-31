%ABSTRACT
% Library function for SDCM. To respect information from perpendicular dimensions (to 
% increase the dynamic range of correlations compared to correlations that are computed 
% only within a signature focus), optionally add flat weights to all genes respectively 
% samples with a total mass of inInfo.searchStrategy.correlation.full2focusWeightRatio
% times the effect size in the twin space. This also stabilizes computation of correlations
% in low-dimensional scenarios and if (nG,nP) is strongly unbalanced (i.e. if 
% min(nG/nP,nP/nG)<<1). (E.g. with nP=2500 samples but only three genes, you should not 
% exclude 2/3 genes from  computation of correlations just because a signature roughly (but 
% not exactly) extends only along one basis gene direction. Instead, focusing should 
% concentrate on the sample space and not the much smaller gene space in this case.)

function WOut = addFlatWeights(W4G,W4P,dimWhereToAddFlatWeights, flatWeights2dynamicWeightsRatio)
  %Initialize:
    assert(sum(size(W4G))>0, 'W4G must be provided; if not yet available, use nan(nG,0)');
    assert(sum(size(W4P))>0, 'W4P must be provided; if not yet available, use nan(0,nP)');
    nG = size(W4G,1);
    nP = size(W4P,2);
  %Get estimated effect sizes (weight sums):
    switch(dimWhereToAddFlatWeights)
      case 1
        effectSizes4genes = sum(abs(W4G),1);
          nCandidates = size(W4G,2);
        if(isempty(W4P))
          effectSizes4samples = nP*ones(1,nCandidates); %assume full-sized effect in sample space, if we know nothing thus far (step 1, only for initial weights).
        else
          effectSizes4samples = sum(abs(W4P),2)'; %support multiple independent candidate vectors along dim 1 of W4P. 
        end
          assert(size(effectSizes4samples,2)==nCandidates, 'W4G and W4P must have the same number of candidate weight vectors (pairs)');
      case 2
        effectSizes4samples = sum(abs(W4P),2);
          nCandidates = size(W4P,1);
        if(isempty(W4G))
          effectSizes4genes = nG*ones(nCandidates,1); %assume a full-sized effect in gene space, if we know nothing thus far (step 1, only for initial weights).
        else
          effectSizes4genes = sum(abs(W4G),1)'; %support multiple independent candidate vectors along dim 2 of W4G. 
        end
          assert(size(effectSizes4genes,1)==nCandidates, 'W4G and W4P must have the same number of candidate weight vectors (pairs)');
    end
  %Define the constant/flat average gene weight :
    switch(dimWhereToAddFlatWeights)
      case 1
        %flatWeights = effectSizes4samples/nG; %effect mass in the twin space distributed uniformly over all dimensions in the current space.
        %effectSizeRatio = effectSizes4genes/nG; %if flatWeight was ==1, scale with this factor so that we add a total weight of effectSizes4genes at the maximum when summed over all nG.
        avgMembership = effectSizes4genes/nG .* (0 + effectSizes4samples/nG); %Note: for numeric 1:1 reproducability, multiply here and leave 0+. (first dividing by nP and multiplying again later loses some digits at precision end... and we did not do this before code cleanup).
      case 2
        %flatWeights = effectSizes4genes/nP; %effect mass in the twin space distributed uniformly over all dimensions in the current space.
        %effectSizeRatio = effectSizes4samples/nP; %if flatWeight was ==1, scale with this factor so that we add a total weight of effectSizes4samples at the maximum when summed over all nP.
        avgMembership = (effectSizes4genes/nP + 0) .* effectSizes4samples/nP; %Note: for numeric 1:1 reproducability, multiply here and leave 0+. (first dividing by nP and multiplying again later loses some digits at precision end... and we did not do this before code cleanup).
    end
  %Add weights to the dimension of the twin space:
    switch(dimWhereToAddFlatWeights)
      case 1
        signs = sign(W4G) + (W4G==0); %prevent zero factors when distributing weights on dimensions of the twin space. In general, the sign of the perpendicular dimension is unknown if its initial effect weight is zero. Without loss of generality, we use a positive weight. Note that correlations etc. use abs(signedW4*) and we want a deterministic solution here.
        %WOut = signs.*bsxfun(@plus, abs(W4G), flatWeights2dynamicWeightsRatio*(effectSizeRatio.*flatWeights));
        WOut = signs.*bsxfun(@plus, abs(W4G), flatWeights2dynamicWeightsRatio*avgMembership);
        WOut = bsxfun(@times, WOut, effectSizes4genes./(eps+sum(abs(WOut),1))); %keep the norm as before adding flat weights; +eps to make 0/0=0 for all zero rows/cols.
        assert(~any(any(isnan(WOut) & ~isnan(W4G))), 'code validation: unexpected NaNs');
      case 2
        signs = sign(W4P) + (W4P==0); %prevent zero factors when distributing weights on dimensions of the twin space. In general, the sign of the perpendicular dimension is unknown if its initial effect weight is zero. Without loss of generality, we use a positive weight. Note that correlations etc. use abs(signedW4*) and we want a deterministic solution here.
        %WOut = signs.*bsxfun(@plus, abs(W4P), flatWeights2dynamicWeightsRatio*(effectSizeRatio.*flatWeights));
        WOut = signs.*bsxfun(@plus, abs(W4P), flatWeights2dynamicWeightsRatio*avgMembership);
        WOut = bsxfun(@times, WOut, effectSizes4samples./(eps+sum(abs(WOut),2))); %keep the norm as before adding flat weights; +eps to make 0/0=0 for all zero rows/cols.
        assert(~any(any(isnan(WOut) & ~isnan(W4P))), 'code validation: unexpected NaNs');
    end
end

