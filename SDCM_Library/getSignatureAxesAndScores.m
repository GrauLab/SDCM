%ABSTRACT
% Library function for SDCM. Collect signature fields from the global 
% results struct (that was filled by precomputeCandidatesInParallel).

function signature = getSignatureAxesAndScores(imj, bOnlyReturnUpToFocusedSignatures)
  %Initialize:
    global eDState;
    signature = eDState.signatureTemplate;
    signature.imj = imj;
    bIsGeneAxisCand = signature.imj>0;
    bIsSampleAxisCand = signature.imj<0;
    if(nargin<2) bOnlyReturnUpToFocusedSignatures = false; end
    if(bIsGeneAxisCand)
      signature.i = +imj;
      eDState.current.precompute.ii = find(eDState.current.precompute.I==signature.imj);
      assert(~isempty(eDState.current.precompute.ii), 'signature for imj=%d was not precomputed, yet', signature.imj);
    end
    if(bIsSampleAxisCand)
      signature.j = -imj;
      eDState.current.precompute.jj = find(eDState.current.precompute.J==-signature.imj);
      assert(~isempty(eDState.current.precompute.jj), 'signature for imj=%d was not precomputed, yet', signature.imj);
    end
  %Get the log2(ratio)s for signature.imj and the twin axes to get a symmetric view on the signature:
    if(bIsGeneAxisCand)
      %Primary axis:
        signature.sampleAxis_origUnits = eDState.current.precompute.sampleAxes4I_origUnits(eDState.current.precompute.ii,:);
        signature.norm4sampleAxis_origUnits = eDState.current.precompute.norm4sampleAxes4I_origUnits(eDState.current.precompute.ii);
      %Twin axis:
        signature.geneAxis_origUnits = eDState.current.precompute.twinGeneAxes4I_origUnits(:,eDState.current.precompute.ii);
        signature.norm4geneAxis_origUnits = eDState.current.precompute.norm4twinGeneAxes4I_origUnits(eDState.current.precompute.ii);
    end
    if(bIsSampleAxisCand)
      %Primary axis:
        signature.geneAxis_origUnits = eDState.current.precompute.geneAxes4J_origUnits(:,eDState.current.precompute.jj);
        signature.norm4geneAxis_origUnits = eDState.current.precompute.norm4geneAxes4J_origUnits(eDState.current.precompute.jj);
      %Twin axis:
        signature.sampleAxis_origUnits = eDState.current.precompute.twinSampleAxes4J_origUnits(eDState.current.precompute.jj,:);
        signature.norm4sampleAxis_origUnits = eDState.current.precompute.norm4twinSampleAxes4J_origUnits(eDState.current.precompute.jj);
    end
    if(bOnlyReturnUpToFocusedSignatures) 
      return; 
    end
  %Correlations of the initial representative candidate to the measured genes respectively samples and the resulting scores:
    %Correlate gene initial representative candidate with all other genes:
      if(bIsGeneAxisCand)
        signature.R4G = eDState.current.precompute.R4Gs4I(:,eDState.current.precompute.ii);
        signature.sampleSizes4R4G = eDState.current.precompute.sampleSizes4R4Gs4I(:,eDState.current.precompute.ii);
        signature.P4R4G = eDState.current.precompute.P4R4Gs4I(:,eDState.current.precompute.ii);
        signature.signedExtendedW4G = eDState.current.precompute.signedExtendedW4Gs4Itwin(:,eDState.current.precompute.ii);
        signature.signedFocusedW4G = eDState.current.precompute.signedFocusedW4Gs4I(:,eDState.current.precompute.ii);
        signature.signedFocusedW4G_withPerpendicularSpace = eDState.current.precompute.signedFocusedW4Gs4I_withPerpendicularSpace(:,eDState.current.precompute.ii);
      end
      if(bIsSampleAxisCand)
        signature.R4G = eDState.current.precompute.R4Gs4Jtwin(:,eDState.current.precompute.jj);
        signature.sampleSizes4R4G = eDState.current.precompute.sampleSizes4R4Gs4Jtwin(:,eDState.current.precompute.jj);
        signature.P4R4G = eDState.current.precompute.P4R4Gs4J(:,eDState.current.precompute.jj);
        signature.signedExtendedW4G = eDState.current.precompute.signedExtendedW4Gs4J(:,eDState.current.precompute.jj);
        signature.signedFocusedW4G = eDState.current.precompute.signedFocusedW4Gs4J(:,eDState.current.precompute.jj);
        signature.signedFocusedW4G_withPerpendicularSpace = eDState.current.precompute.signedFocusedW4Gs4J_withPerpendicularSpace(:,eDState.current.precompute.jj);
      end
    %Correlate sample initial representative candidate with all other samples:
      if(bIsGeneAxisCand)
        signature.R4P = eDState.current.precompute.R4Ps4Itwin(eDState.current.precompute.ii,:);
        signature.sampleSizes4R4P = eDState.current.precompute.sampleSizes4R4Ps4Itwin(eDState.current.precompute.ii,:);
        signature.P4R4P = eDState.current.precompute.P4R4Ps4I(eDState.current.precompute.ii,:);
        signature.signedExtendedW4P = eDState.current.precompute.signedExtendedW4Ps4I(eDState.current.precompute.ii,:);
        signature.signedFocusedW4P = eDState.current.precompute.signedFocusedW4Ps4I(eDState.current.precompute.ii,:);
        signature.signedFocusedW4P_withPerpendicularSpace = eDState.current.precompute.signedFocusedW4Ps4I_withPerpendicularSpace(eDState.current.precompute.ii,:);
      end
      if(bIsSampleAxisCand)
        signature.R4P = eDState.current.precompute.R4Ps4J(eDState.current.precompute.jj,:);
        signature.sampleSizes4R4P = eDState.current.precompute.sampleSizes4R4Ps4J(eDState.current.precompute.jj,:);
        signature.P4R4P = eDState.current.precompute.P4R4Ps4J(eDState.current.precompute.jj,:);
        signature.signedExtendedW4P = eDState.current.precompute.signedExtendedW4Ps4Jtwin(eDState.current.precompute.jj,:);
        signature.signedFocusedW4P = eDState.current.precompute.signedFocusedW4Ps4J(eDState.current.precompute.jj,:);
        signature.signedFocusedW4P_withPerpendicularSpace = eDState.current.precompute.signedFocusedW4Ps4J_withPerpendicularSpace(eDState.current.precompute.jj,:);
      end
    %Signature size estimates and statistics:
      if(bIsGeneAxisCand)
        signature.signatureSizeByCorrSum2D = eDState.current.precompute.signatureSizeByCorrSum2D4I(eDState.current.precompute.ii);
        signature.signatureSizeByCorrSum4G = eDState.current.precompute.signatureSizeByCorrSum4G4I(eDState.current.precompute.ii);
        signature.signatureSizeByCorrSum4P = eDState.current.precompute.signatureSizeByCorrSum4P4I(eDState.current.precompute.ii);
        signature.signatureAbsMean2D = eDState.current.precompute.signatureAbsMean2D4I(eDState.current.precompute.ii);
        signature.signatureCorrInExtendedFocus = eDState.current.precompute.signatureCorrInExtendedFocus4I(eDState.current.precompute.ii);
        signature.sampleSize4signatureAbsMean2D = eDState.current.precompute.sampleSize4signatureAbsMean2D4I(eDState.current.precompute.ii);
        signature.signatureAbsSD2D = eDState.current.precompute.signatureAbsSD2D4I(eDState.current.precompute.ii);
        signature.log10_p = eDState.current.precompute.log10_p4I(eDState.current.precompute.ii);
        signature.log10_p4Correlations = eDState.current.precompute.log10_p4Correlations4I(eDState.current.precompute.ii);
        signature.log10_p4SignalStrength = eDState.current.precompute.log10_p4SignalStrength4I(eDState.current.precompute.ii);
      end
      if(bIsSampleAxisCand)
        signature.signatureSizeByCorrSum2D = eDState.current.precompute.signatureSizeByCorrSum2D4J(eDState.current.precompute.jj);
        signature.signatureSizeByCorrSum4G = eDState.current.precompute.signatureSizeByCorrSum4G4J(eDState.current.precompute.jj);
        signature.signatureSizeByCorrSum4P = eDState.current.precompute.signatureSizeByCorrSum4P4J(eDState.current.precompute.jj);
        signature.signatureAbsMean2D = eDState.current.precompute.signatureAbsMean2D4J(eDState.current.precompute.jj);
        signature.signatureCorrInExtendedFocus = eDState.current.precompute.signatureCorrInExtendedFocus4J(eDState.current.precompute.jj);
        signature.sampleSize4signatureAbsMean2D = eDState.current.precompute.sampleSize4signatureAbsMean2D4J(eDState.current.precompute.jj);
        signature.signatureAbsSD2D = eDState.current.precompute.signatureAbsSD2D4J(eDState.current.precompute.jj);
        signature.log10_p = eDState.current.precompute.log10_p4J(eDState.current.precompute.jj);
        signature.log10_p4Correlations = eDState.current.precompute.log10_p4Correlations4J(eDState.current.precompute.jj);
        signature.log10_p4SignalStrength = eDState.current.precompute.log10_p4SignalStrength4J(eDState.current.precompute.jj);
      end
    assert(isempty(setdiff(fieldnames(signature),fieldnames(eDState.signatureTemplate))), 'code validation: the signature template does not include all actually used signature fields');
end

