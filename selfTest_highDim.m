%ABSTRACT
%  Demo and self test of SDCM for a high-dimensional test signal
%  simulating thirteen distinct interactions.
%SYNTAX
%  Just run this script (press F5).

%% Initialize:
  %Warn on missing dependencies:
    if(true)
      installedToolboxesInfo = ver;
      csRequiredToolboxes = {
        {'Statistics Toolbox','Statistics and Machine Learning Toolbox'} %'Statistics Toolbox' is the older name of this toolbox
        {'Bioinformatics Toolbox'}
        {'Parallel Computing Toolbox'}
      };
      for i=1:length(csRequiredToolboxes)
        if(~any(ismember(csRequiredToolboxes{i}, {installedToolboxesInfo.Name})))
          if(isscalar(csRequiredToolboxes{i}))
            sMissingToolbox = csRequiredToolboxes{i}{1};
          else
            sMissingToolbox = ['with a name in {',strtrim(evalc('disp(csRequiredToolboxes{i}(:))')),'}'];
          end
          warning('Dependencies check: MISSING TOOLBOX %s\n(errors will occur, once functions of this toolbox are invoked by SDCM)', sMissingToolbox);
        end
      end
    end
  %Check access to library functions:
    if(exist('SDCM_printStatus','file')~=2)
      sLibrarySubdir = 'SDCM_Library';
      sPathToSelf = fileparts(mfilename('fullpath'));
      if(exist(fullfile(sPathToSelf,sLibrarySubdir,'tools','SDCM_printStatus.m'),'file')>0)
        warning('Now adding SDCM library functions in [%s] to the Matlab path.', fullfile(sPathToSelf,sLibrarySubdir));
        addpath(genpath(fullfile(sPathToSelf,sLibrarySubdir)));
        assert(exist('SDCM_printStatus','file')==2, 'Cannot find SDCM library functions such as SDCM_printStatus.m even after adding %s to the Matlab path. Please check your installation and add the library functions manually to your Matlab path.', fullfile(sPathToSelf,sLibrarySubdir));
      else
        error('Cannot find SDCM library functions subfolder at %s. Please add the library functions manually to your Matlab path. Pressing F5 and then changing the directory to this self test file should also be sufficient.', sPathToSelf);
      end
    end

%% Simulation of a high-dimensional test signal containing several superposed signatures with random gene and sample orders:
  %Simulate data:
    fprintf('-> simulating data ...\n');
    simulated = struct();
    simulated.nG = 1000;
    simulated.nP = 100;
    simulated.rowIDs = cellfun(@(n)sprintf('i%d',n),num2cell((1:simulated.nG)'),'UniformOutput',false);
    simulated.colIDs = cellfun(@(n)sprintf('j%d',n),num2cell((1:simulated.nP)),'UniformOutput',false);
    simulated.signatures = struct();
    %simulate measurement noise with sigma=1/2:
      simulated.noiseSigma = 1/2;
      %simulated.noiseSigma = 1;
      simulated.simNoise = simulated.noiseSigma*randn(simulated.nG,simulated.nP);
    %simulate missing values:
      ratioOfMissingValues = 0; %keep zero for the comprison plot (PCA cannot handle missing values...)
      %ratioOfMissingValues = 0.4; %0.1; 
      if(ratioOfMissingValues>0)
        linI = randperm(numel(simulated.simNoise));
        simulated.simNoise(linI(1:ceil(numel(simulated.simNoise)*ratioOfMissingValues))) = NaN;
      end
    %Simulate all true signature signals configured for the current simulationScenario:
      %Select the test scenario to simulate:
        sTestMnemonic = 'superposition test with 7x pattern #3'; %Hint: to push detection limits try increasing ratioOfMissingValues to 40% and double .noiseSigma above.
        %sTestMnemonic = 'versatility test with 19 signatures'; %simulates 19 superposed signatures of 7 distinct patterns as described in the paper for the smaller/easier versatility test.
          %<-Note: Signatures of the weak pattern #6 are sometimes not detected in the versa19 test. Due to the relatively high superposition depth, their very weak signal does not sufficiently stand out of the estimated noise level. Read the "best scores seen" part at the end of the detection log for details. To detect these weak patterns anyway, you can try reducing qualification thresholds (e.g. reduce configuration.searchStrategy.qualification.minCorrInExtendedFocus=0.33 and only demand configuration.searchStrategy.qualification.alpha4signalStrength=1e-3). However, reducing false negatives in this way may generally invite false positives.
        fprintf('    -> selected test scenario: "%s"...\n', sTestMnemonic);
      fcnMeanW = @weightedMean_supportsNaNs; %dissectSignal('getLibraryFunction','fcnMeanW_supportsNans');
      fcnRandPerm = @(n)randperm(n);
      IEligibleGenes = 1:simulated.nG;
      k=0;
      while(true) k=k+1;
        %Declare signature to simulate:
          switch(sTestMnemonic)
            case 'versatility test with 7 signatures'
              switch(k)
                case 1
                  simulated.signatures(k).mnemonic = 'gradual from strong to nonexisting in genes, binary strong in samples\nall genes affected (1:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'technical signatures like different labeling protocols can cause very broad and strong overlay signatures (easy to detect, important to remove)';
                  simulated.signatures(k).ratioOfAffectedGenes = 1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1.5;
                  simulated.signatures(k).minL2R = 0;
                  simulated.signatures(k).fcnNormalizedSignaturePos2Folding = @(posG,posP)...
                      simulated.signatures(k).maxL2R ...
                    .*sign(posG-simulated.signatures(k).ratioOfTGZeroTransition+eps)...
                    .*sign(posP-simulated.signatures(k).ratioOfTPZeroTransition+eps) ...
                    .*abs(ilv(...
                        posG-simulated.signatures(k).ratioOfTGZeroTransition... %just gradual over genes, but sharp between patient partitions.
                       ,@(S) simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R...
                            +(1-simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R) * S/max(max(abs(S)))...
                      )) ...
                  ;
                case {2}
                  simulated.signatures(k).mnemonic = 'binary and strong in both genes and samples\n1%% of all genes affected (7:3), all samples affected (1:1)';
                  simulated.signatures(k).example = 'gender clusters typically show strong binary signals for few measured genes';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.01;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 7/(7+3);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 2;
                case {3}
                  simulated.signatures(k).mnemonic = 'gradual, strong to weak\n20%% of all genes affected (3:1), 50%% of all samples affected (1:2)';
                  simulated.signatures(k).example = 'broad biological signature that only exists in a subclass of samples, in some stronger in others weaker';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.2;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 3/(3+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+2);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1/2;
                case {4}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n10%% of all genes affected (0:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'medium biological signature';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {5}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n2%% of all genes affected (1:0), 50%% of all samples affected (1:1)';
                  simulated.signatures(k).example = 'medium narrow biological signature that can only increase expressions above normal';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.02;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1;
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {6}
                  simulated.signatures(k).mnemonic = 'gradual, weak to nonexisting\n10%% of all genes affected (1:1), all samples affected (1:1)';
                  simulated.signatures(k).example = 'very weak biological signature (maximum signal is 1\\sigma of the simulated normal noise)';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1/2;
                  simulated.signatures(k).minL2R = 0;
                case {7}
                  simulated.signatures(k).mnemonic = 'gradual, strong to medium\n25%% of all genes affected (0:1), only 5%% of all samples affected (0:1)';
                  simulated.signatures(k).example = 'strong biological signature in only few samples, typical for immunities';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.25;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.05;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 0;
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1;
                otherwise
                  break;
              end
              fprintf('    <- simulated pattern %d/%d\n', k, 7);
            case 'versatility test with 13 signatures'
              switch(k)
                case 1
                  simulated.signatures(k).mnemonic = 'gradual from strong to nonexisting in genes, binary strong in samples\nall genes affected (1:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'technical signatures like different labeling protocols can cause very broad and strong overlay signatures (easy to detect, important to remove)';
                  simulated.signatures(k).ratioOfAffectedGenes = 1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1.5;
                  simulated.signatures(k).minL2R = 0;
                  simulated.signatures(k).fcnNormalizedSignaturePos2Folding = @(posG,posP)...
                      simulated.signatures(k).maxL2R ...
                    .*sign(posG-simulated.signatures(k).ratioOfTGZeroTransition+eps)...
                    .*sign(posP-simulated.signatures(k).ratioOfTPZeroTransition+eps) ...
                    .*abs(ilv(...
                        posG-simulated.signatures(k).ratioOfTGZeroTransition... %just gradual over genes, but sharp between patient partitions.
                       ,@(S) simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R...
                            +(1-simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R) * S/max(max(abs(S)))...
                      )) ...
                  ;
                case {2,8,11}
                  simulated.signatures(k).mnemonic = 'binary and strong in both genes and samples\n1%% of all genes affected (7:3), all samples affected (1:1)';
                  simulated.signatures(k).example = 'gender clusters typically show strong binary signals for few measured genes';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.01;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 7/(7+3);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 2;
                case {3,9,12}
                  simulated.signatures(k).mnemonic = 'gradual, strong to weak\n20%% of all genes affected (3:1), 50%% of all samples affected (1:2)';
                  simulated.signatures(k).example = 'broad biological signature that only exists in a subclass of samples, in some stronger in others weaker';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.2;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 3/(3+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+2);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1/2;
                case {4,10,13}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n10%% of all genes affected (0:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'medium biological signature';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {5}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n2%% of all genes affected (1:0), 50%% of all samples affected (1:1)';
                  simulated.signatures(k).example = 'medium narrow biological signature that can only increase expressions above normal';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.02;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1;
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {6}
                  simulated.signatures(k).mnemonic = 'gradual, weak to nonexisting\n10%% of all genes affected (1:1), all samples affected (1:1)';
                  simulated.signatures(k).example = 'very weak biological signature (maximum signal is 1\\sigma of the simulated normal noise)';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1/2;
                  simulated.signatures(k).minL2R = 0;
                case {7}
                  simulated.signatures(k).mnemonic = 'gradual, strong to medium\n25%% of all genes affected (0:1), only 5%% of all samples affected (0:1)';
                  simulated.signatures(k).example = 'strong biological signature in only few samples, typical for immunities';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.25;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.05;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 0;
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1;
                otherwise
                  break;
              end
              fprintf('    <- simulated pattern %d/%d\n', k, 13);
            case 'versatility test with 19 signatures'
              switch(k)
                case 1
                  simulated.signatures(k).mnemonic = 'gradual from strong to nonexisting in genes, binary strong in samples\nall genes affected (1:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'technical signatures like different labeling protocols can cause very broad and strong overlay signatures (easy to detect, important to remove)';
                  simulated.signatures(k).ratioOfAffectedGenes = 1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1.5;
                  simulated.signatures(k).minL2R = 0;
                  simulated.signatures(k).fcnNormalizedSignaturePos2Folding = @(posG,posP)...
                      simulated.signatures(k).maxL2R ...
                    .*sign(posG-simulated.signatures(k).ratioOfTGZeroTransition+eps)...
                    .*sign(posP-simulated.signatures(k).ratioOfTPZeroTransition+eps) ...
                    .*abs(ilv(...
                        posG-simulated.signatures(k).ratioOfTGZeroTransition... %just gradual over genes, but sharp between patient partitions.
                       ,@(S) simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R...
                            +(1-simulated.signatures(k).minL2R/simulated.signatures(k).maxL2R) * S/max(max(abs(S)))...
                      )) ...
                  ;
                case {2,8,14}
                  simulated.signatures(k).mnemonic = 'binary and strong in both genes and samples\n1%% of all genes affected (7:3), all samples affected (1:1)';
                  simulated.signatures(k).example = 'gender clusters typically show strong binary signals for few measured genes';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.01;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 7/(7+3);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 2;
                case {3,9,15}
                  simulated.signatures(k).mnemonic = 'gradual, strong to weak\n20%% of all genes affected (3:1), 50%% of all samples affected (1:2)';
                  simulated.signatures(k).example = 'broad biological signature that only exists in a subclass of samples, in some stronger in others weaker';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.2;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 3/(3+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+2);
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1/2;
                case {4,10,16}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n10%% of all genes affected (0:1), all samples affected (2:1)';
                  simulated.signatures(k).example = 'medium biological signature';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 2/(2+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {5,11,17}
                  simulated.signatures(k).mnemonic = 'gradual, medium to weak\n2%% of all genes affected (1:0), 50%% of all samples affected (1:1)';
                  simulated.signatures(k).example = 'medium narrow biological signature that can only increase expressions above normal';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.02;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1;
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1;
                  simulated.signatures(k).minL2R = 1/2;
                case {6,12,18}
                  simulated.signatures(k).mnemonic = 'gradual, weak to nonexisting\n10%% of all genes affected (1:1), all samples affected (1:1)';
                  simulated.signatures(k).example = 'very weak biological signature (maximum signal is 1\\sigma of the simulated normal noise)';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.1;
                  simulated.signatures(k).ratioOfAffectedSamples = 1;
                  simulated.signatures(k).ratioOfTGZeroTransition = 1/(1+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+1);
                  simulated.signatures(k).maxL2R = 1/2;
                  simulated.signatures(k).minL2R = 0;
                case {7,13,19}
                  simulated.signatures(k).mnemonic = 'gradual, strong to medium\n25%% of all genes affected (0:1), only 5%% of all samples affected (0:1)';
                  simulated.signatures(k).example = 'strong biological signature in only few samples, typical for immunities';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.25;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.05;
                  simulated.signatures(k).ratioOfTGZeroTransition = 0;
                  simulated.signatures(k).ratioOfTPZeroTransition = 0;
                  simulated.signatures(k).maxL2R = 2;
                  simulated.signatures(k).minL2R = 1;
                otherwise
                  break;
              end
              fprintf('    <- simulated pattern %d/%d\n', k, 19);
            case 'superposition test with 7x pattern #3'
              switch(k)
                case num2cell(1:7)
                %case num2cell(1:16)
                %case num2cell(1:20)
                  simulated.signatures(k).mnemonic = 'gradual, strong to weak\n20%% of all genes affected (3:1), 50%% of all samples affected (1:2)';
                  simulated.signatures(k).example = 'broad biological signature that only exists in a subclass of samples, in some stronger in others weaker';
                  simulated.signatures(k).ratioOfAffectedGenes = 0.2;
                  simulated.signatures(k).ratioOfAffectedSamples = 0.5;
                  simulated.signatures(k).ratioOfTGZeroTransition = 3/(3+1);
                  simulated.signatures(k).ratioOfTPZeroTransition = 1/(1+2);
                  simulated.signatures(k).maxL2R = 2 *(-1)^k; %times (-1)^k to not accumulate a global offset effect by superposition (the signature template is asymetric).
                  simulated.signatures(k).minL2R = 1/2 *(-1)^k; %times (-1)^k to not accumulate a global offset effect by superposition (the signature template is asymetric).
                otherwise
                  break;
              end
              fprintf('    <- simulated pattern %d/%d\n', k, 7);
          end
        %Get the affected genes and samples in a randomly simulated signature:
          if(true)
            %randomly select affected genes and samples:
              RI = IEligibleGenes(fcnRandPerm(length(IEligibleGenes))); %Support restricting signatures simulation to a subspace (leaving (1-signatureGeneSizesScaleFactor)*simulated.nG noise genes).
              simulated.signatures(k).SI = RI(1:round(length(RI)*simulated.signatures(k).ratioOfAffectedGenes))';
              RJ = fcnRandPerm(simulated.nP);
              simulated.signatures(k).SJ = RJ(1:round(length(RJ)*simulated.signatures(k).ratioOfAffectedSamples));
            %order the affected samples on a eqd axis:
              posG = linspace(0,1,length(simulated.signatures(k).SI))';
              posP = linspace(0,1,length(simulated.signatures(k).SJ));
            %Save affected samples/genes for later:
              simulated.signatures(k).affectedSamples = zeros(1,simulated.nP);
                simulated.signatures(k).affectedSamples(simulated.signatures(k).SJ) = sign(eps+posP-simulated.signatures(k).ratioOfTPZeroTransition);
                if(sum(simulated.signatures(k).affectedSamples~=0)<2)
                  warning('less than 2 samples affected by simulated signature %d; it principally cannot be detected',k);
                end
              simulated.signatures(k).affectedGenes = zeros(simulated.nG,1);
                simulated.signatures(k).affectedGenes(simulated.signatures(k).SI) = sign(eps+posG-simulated.signatures(k).ratioOfTGZeroTransition);
                if(sum(simulated.signatures(k).affectedGenes~=0)<2)
                  warning('less than 2 genes affected by simulated signature %d; it principally cannot be detected',k);
                end
          end
        %Simulate foldings for current signature:
          if(true)
            simulated.signatures(k).Pk = zeros(simulated.nG,simulated.nP);
            %Set common general fcn if no custom has been provided:
              if(~isfield(simulated.signatures(k),'fcnNormalizedSignaturePos2Folding') || isempty(simulated.signatures(k).fcnNormalizedSignaturePos2Folding))
                %simulated.signatures(k).exponent = 1; %0.1;
                simulated.signatures(k).fcnNormalizedSignaturePos2Folding = @(posG,posP)...
                    simulated.signatures(k).maxL2R ...
                  .*sign(posG-simulated.signatures(k).ratioOfTGZeroTransition+eps)...
                  .*sign(posP-simulated.signatures(k).ratioOfTPZeroTransition+eps) ...
                  .*abs(ilv(...
                      posG-simulated.signatures(k).ratioOfTGZeroTransition...
                     ,posP-simulated.signatures(k).ratioOfTPZeroTransition...
                     ,@(SG,SP)sqrt(abs(SG.*SP))...
                     ,@(S)  simulated.signatures(k).minL2R / simulated.signatures(k).maxL2R ...
                          + (1 - simulated.signatures(k).minL2R / simulated.signatures(k).maxL2R) * S/max(max(abs(S))) ...
                     ...,@(L2Rs) sign(L2Rs).*abs(L2Rs).^simulated.signatures(k).exponent ...
                    )) ...
                ;
              end
            %Compute the configured bimonotonic signal and assign to Pk in reference gene and sample orders:
              simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ) = simulated.signatures(k).fcnNormalizedSignaturePos2Folding(...
                repmat(posG,1,length(simulated.signatures(k).SJ)) ...
               ,repmat(posP,length(simulated.signatures(k).SI),1) ...
              );
              %<-Assert monotonicity of declared function:
                simulated.signatures(k).BPAllIncreasing = all(diff(simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ),1,1)>=0,1);
                simulated.signatures(k).BPAllDecreasing = all(diff(simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ),1,1)<=0,1);
                  assert(all(simulated.signatures(k).BPAllIncreasing|simulated.signatures(k).BPAllDecreasing), 'some samples are not monotonic over gene axis in signature');
                simulated.signatures(k).BGAllIncreasing = all(diff(simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ),1,2)>=0,2);
                simulated.signatures(k).BGAllDecreasing = all(diff(simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ),1,2)<=0,2);
                  assert(all(simulated.signatures(k).BGAllIncreasing|simulated.signatures(k).BGAllDecreasing), 'some genes are not monotonic over sample axis in signature');
          end
      end
      %Compute gene axes and correlations (R4G, R4P) for simulated signatures (for later comparison):
        fprintf('    -> computing gene/sample axes for simulated signatures...\n');
        for k=1:length(simulated.signatures)
          %Infer axes from the exact simulated signature eigensignal (without noise), simply by using the exactly known affected genes and samples for focusing:
            L2Rs = simulated.signatures(k).Pk; %use exact/true sampleAcnhors without noise for later comparison of detected results with simulated signatures.
            simulated.signatures(k).sampleAxis_origUnits = fcnMeanW(...
              L2Rs ...
             ,simulated.signatures(k).affectedGenes... %we know the focus exactly for simulated patterns; so use it as initial weights.
            ,1);
            simulated.signatures(k).geneAxis_origUnits = fcnMeanW(...
              L2Rs ...
             ,simulated.signatures(k).affectedSamples... %we know the focus exactly for simulated patterns; so use it as initial weights.
            ,2);
        end
      %Add rowIDs and colIDs:
        fprintf('    -> defining row and column IDs ...\n');
        for k=1:length(simulated.signatures)
          B = simulated.signatures(k).affectedGenes==+1;
            simulated.rowIDs(B) = cellfun(@(s)sprintf('#%d%s,%s',k,'+',s),simulated.rowIDs(B),'UniformOutput',false);
          B = simulated.signatures(k).affectedGenes==-1;
            simulated.rowIDs(B) = cellfun(@(s)sprintf('#%d%s,%s',k,'-',s),simulated.rowIDs(B),'UniformOutput',false);
          B = simulated.signatures(k).affectedSamples==+1;
            simulated.colIDs(B) = cellfun(@(s)sprintf('#%d%s,%s',k,'+',s),simulated.colIDs(B),'UniformOutput',false);
          B = simulated.signatures(k).affectedSamples==-1;
            simulated.colIDs(B) = cellfun(@(s)sprintf('#%d%s,%s',k,'-',s),simulated.colIDs(B),'UniformOutput',false);
        end
      %superpose all signatures:
        simulated.simL2Rs = zeros(simulated.nG,simulated.nP);
        for k=1:length(simulated.signatures)
          simulated.simL2Rs = simulated.simL2Rs + simulated.signatures(k).Pk;
        end
  %Plot simulated signatures:
    %initialize colormaps for plots:
      if(true) %CM for heatmaps
        colorResolution=254;
        ramp = ((1:colorResolution/2)/(colorResolution/2))';  
        CM = flipud([[flipud(ramp) zeros(size(ramp)) zeros(size(ramp))]; [0 0 0]; [zeros(size(ramp)) zeros(size(ramp)), ramp]]); %Nature Methods does not want red/green.
        CM = bsxfun(@times,sum(CM,2).^(1/3),(1-bsxfun(@times,1-sum(CM,2)/3,exp(-3*CM))).^1); %mehr dynamic range.
      end
      if(true) %CM4Rs
        colorResolution=250;
        ramp = ((1:colorResolution)/(colorResolution))';  
        CM4Rs = [ramp,ramp,ramp];
        CM4Rs(end/2:end,2:3) = 1-CM4Rs(end/2:end,2:3);
      end
    %Plot each simulated signature in its ordered/bimotonic pattern:
      if(true)
        pixelFactor = 250; 
        bShowNaNsInSimulatedSignatures = ratioOfMissingValues>0;
        for k=1:length(simulated.signatures)
          f=figure('Position',[73 1 505 997],'Color','w');
            CM_supportsNaNs = [[.5 .5 .5];CM]; cLim=[-2,2]; hold on;
              colormap(CM_supportsNaNs); caxis(cLim); 
            data2Plot = simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ);
              if(bShowNaNsInSimulatedSignatures)
                data2Plot(isnan(simulated.simNoise(simulated.signatures(k).SI,simulated.signatures(k).SJ))) = NaN; %add simulated NaNs.
              end
              lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(CM,1);
                data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
              imagesc(data2Plot);
              daspect([1 1 1]);
              axis tight;
            title(sprintf('simulated signature %d\nmax signal = %0.2f', k,max(max(abs(simulated.signatures(k).Pk(simulated.signatures(k).SI,simulated.signatures(k).SJ))))),'FontSize',8);
              set(gca,'XTick',[],'YTick',[]);
              xlabel(sprintf('%d+%d samples',sum(simulated.signatures(k).affectedSamples==+1),sum(simulated.signatures(k).affectedSamples==-1)));
              ylabel(sprintf('%d+%d genes',sum(simulated.signatures(k).affectedGenes==+1),sum(simulated.signatures(k).affectedGenes==-1)));
            set(gca,'Unit','Pixels','Position',[50, 50, sqrt(pixelFactor*length(simulated.signatures(k).SJ)), sqrt(pixelFactor*length(simulated.signatures(k).SI))]);
          drawnow;
        end
      end
    %Plot the superposition of all simualted signatures (input signal):
      if(true)
        f=figure('Position',[73 1 505 997],'Color','w');
          CM_supportsNaNs = [[.5 .5 .5];CM]; cLim=[-2,2];
            colormap(CM_supportsNaNs); caxis(cLim); hold on;
          data2Plot = simulated.simL2Rs+simulated.simNoise;
            lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(CM,1);
              data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
            imagesc(data2Plot);
            daspect([1 1 1]);
            axis tight;
          title(sprintf('sum of all simulated signatures\nplus simulated normal noise (\\sigma = %0.2f log_2(ratio)s) and %1.0f%% missing values',simulated.noiseSigma,100*ratioOfMissingValues),'FontSize',8);
            set(gca,'XTick',[],'YTick',[]);
            xlabel(sprintf('%d samples',simulated.nP));
            ylabel(sprintf('%d genes',simulated.nG));
          cbh=colorbar;
            if(true) %remove the gray color from the color bar:
              hImg = findobj(cbh,'Type','Image');
              if(~isempty(hImg)) %colormap image no longer accessible since Matlab's new graphics engine.
                CData = get(hImg,'CData');
                if(CData(end)==size(CM_supportsNaNs,1))
                  CData(1)=2;
                end
                set(findobj(cbh,'Type','Image'),'CData',CData);
              end
            end
      end
      
%% Dissect the signal:
  configuration = SDCM_defaultConfig(simulated.nG, simulated.nP);
    configuration.export.nStatusOutputLevel = 1; %disable most status messages.
  %Plot config:
    configuration.plots.bVisible = true; %show plots for user.          
    configuration.reference.rowIDs = simulated.rowIDs;
    configuration.reference.colIDs = simulated.colIDs;
    configuration.plots.descriptions.sTitlePrefix = 'detected signature ';
    configuration.plots.descriptions.sSamplesMnemonic = 'samples';
    configuration.plots.descriptions.sGenesMnemonic = 'genes';
    configuration.plots.descriptions.sSamplesLabel = sprintf('simulated samples (of %d)', simulated.nP);
    configuration.plots.descriptions.sGenesLabel = sprintf('simulated genes (of %d)', simulated.nG);
    configuration.plots.signatureDefinitions.cLim4L2Rs = [-2.5,2.5]; %use a common non-adaptive color axis for all signatures in the versatility test to see signal strength differences.
    configuration.plots.focusConvergence.bEnabled = true;
  %config modification example: increase detection sensitivity by reducing default qualification thresholds (useful for detection of the weak signal pattern #6 even in context of the high superposition-depth 'versatility test with 19 signatures'; however, this may invite some false positive signatures):
    %configuration.searchStrategy.qualification.alpha4signalStrength = 1e-3; %default=1e-5
    %configuration.searchStrategy.qualification.minCorrInExtendedFocus = 0.33; %default=0.4
  %Start parallelization: 
    nCores = 4; %set to the number of cores your machine has.
    if(isempty(gcp)) parpool(nCores); end
  %Call dissection:
    inputSignal = simulated.simL2Rs+simulated.simNoise;
    detection = SDCM(inputSignal, configuration);
  %For a comparison, compute principal components:
    detectionPCA = struct();
      [detectionPCA.sampleAxes,~,detectionPCA.variancesOfsampleAxes] = pca(inputSignal'); %,'Economy',false);
      [detectionPCA.geneAxes,~,detectionPCA.variancesOfgeneAxes] = pca(inputSignal); %,'Economy',false);
        detectionPCA.geneAxes = detectionPCA.geneAxes';
    if(any(isnan(inputSignal)))
      warning('PCA does not support missing values; for comparison with PCA set ratioOfMissingValues=0 above');
    end

%% Correlate detection results with simulated signature axes:
  for sComparisonTarget = {'methodValidation','PCA'}; sComparisonTarget=sComparisonTarget{1};
    %Correlate simulated and detected signature axes:
      for k=1:length(simulated.signatures)
        %initialize:
          switch(sComparisonTarget) 
            case 'methodValidation'
              nDetected = length(detection.signatures);
              simulated.signatures(k).corr2detectedGeneAxis = nan(nDetected,1);
              simulated.signatures(k).corr2detectedSampleAxis = nan(nDetected,1);
            case 'PCA'
              %nDetected = size(detectionPCA.sampleAxes,2); 
              nDetected = length(simulated.signatures); %PCA has no concept of false positives and returns as many PCs as the space has dimensions; we only plot the strongest PCs by variance via providing PCA the info of the actual number of simulated signatures here; use the previous line to plot all PCs.
              simulated.signatures(k).corr2detectedGeneSpacePC = nan(nDetected,1);
              simulated.signatures(k).corr2detectedSampleSpacePC = nan(nDetected,1);
          end
        %Correlate with all detected axes:
          for l=1:nDetected
            %Gene axes comparisons:
              if(true)
                %Get simulated and detected gene axes to compare:
                  sourceSignature4G = simulated.signatures(k).geneAxis_origUnits;
                  switch(sComparisonTarget)
                    case 'methodValidation'
                      targetSignature4G = detection.signatures(l).step2_finalSignatureAxes.geneAxis_origUnits;
                    case 'PCA'
                      targetSignature4G = nan(size(simulated.signatures(k).geneAxis_origUnits));
                      if(size(detectionPCA.sampleAxes,2)>=l) %can be empty for NaN data
                        targetSignature4G = detectionPCA.sampleAxes(:,l);
                      end
                  end
                %Define weights: As PCA does not provide dimension weights, use a common weights formula based only on the signal for a fair comparison (instead of using the SDCM signature focus):
                  W4Gsim = ilv(simulated.signatures(k).geneAxis_origUnits, @(X)min(1, abs(X)/max(abs(X))/0.5)); %relative signal, 50% of max(abs(signal)) already maps to 100% weight.
                  W4Gdet = ilv(targetSignature4G, @(X)min(1, abs(X)/max(abs(X))/0.5)); %relative signal, 50% of max(abs(signal)) already maps to 100% weight.
                  W4G = max(W4Gsim, W4Gdet); %compare all dims that are either strong in the simulated or in the detected effect.
                %Correlate:
                  signatureCorrs4G = uncenteredWeightedCorrelation(...
                     sourceSignature4G, targetSignature4G, 1 ...
                    ,W4G ...
                  );
                %save:
                  switch(sComparisonTarget) 
                    case 'methodValidation'
                      simulated.signatures(k).corr2detectedGeneAxis(l) = abs(signatureCorrs4G);
                      detection.signatures(l).corr2simulatedGeneAxis(k) = abs(signatureCorrs4G);
                    case 'PCA'
                      simulated.signatures(k).corr2detectedGeneSpacePC(l) = abs(signatureCorrs4G);
                      detectionPCA.PCs(l).corr2simulatedGeneSpacePC(k) = abs(signatureCorrs4G);
                  end
              end
            %Sample axes comparisons:
              if(true)
                %Get simulated and detected sample axes to compare:
                  sourceSignature4P = simulated.signatures(k).sampleAxis_origUnits;
                  switch(sComparisonTarget)
                    case 'methodValidation'
                      targetSignature4P = detection.signatures(l).step2_finalSignatureAxes.sampleAxis_origUnits;
                    case 'PCA'
                      targetSignature4P = nan(size(simulated.signatures(k).sampleAxis_origUnits));
                      if(size(detectionPCA.geneAxes,1)>=l) %can be empty for NaN data
                        targetSignature4P = detectionPCA.geneAxes(l,:);
                      end
                  end
                %Define weights: As PCA does not provide dimension weights, use a common weights formula based only on signal strengths for a fair comparison (instead of using the SDCM signature focus):
                  W4Psim = min(1, abs(simulated.signatures(k).sampleAxis_origUnits)/max(abs(simulated.signatures(k).sampleAxis_origUnits))/0.5); %relative signal, 50% of max(abs(signal)) already maps to 100% weight.
                  W4Pdet = min(1, abs(targetSignature4P)/max(abs(targetSignature4P))/0.5); %relative signal, 50% of max(abs(signal)) already maps to 100% weight.
                  W4P = max(W4Psim, W4Pdet); %compare all dims that are either strong in the simulated or in the detected effect.
                %Correlate:
                  signaturesCorr4P = uncenteredWeightedCorrelation(...
                     sourceSignature4P, targetSignature4P, 2 ...
                    ,W4P ...
                  );
                %save:
                  switch(sComparisonTarget) 
                    case 'methodValidation'
                      simulated.signatures(k).corr2detectedSampleAxis(l) = abs(signaturesCorr4P);
                      detection.signatures(l).corr2simulatedSampleAxis(k) = abs(signaturesCorr4P);
                    case 'PCA'
                      simulated.signatures(k).corr2detectedSampleSpacePC(l) = abs(signaturesCorr4P);
                      detectionPCA.PCs(l).corr2simulatedSampleSpacePC(k) = abs(signaturesCorr4P);
                  end
              end
          end
      end
    %Collect all simulated<->detected axes correlations:
      switch(sComparisonTarget)
        case 'methodValidation'
          Rs4detectedAndSimulatedGeneAxes = zeros(nDetected, length(simulated.signatures));
          Rs4detectedAndSimulatedSampleAxes = zeros(nDetected, length(simulated.signatures));
          for l=1:nDetected
            Rs4detectedAndSimulatedGeneAxes(l,:) = detection.signatures(l).corr2simulatedGeneAxis;
            Rs4detectedAndSimulatedSampleAxes(l,:) = detection.signatures(l).corr2simulatedSampleAxis;
          end
        case 'PCA'
          Rs4detectedAndSimulatedGeneAxes = zeros(nDetected, length(simulated.signatures));
          Rs4detectedAndSimulatedSampleAxes = zeros(nDetected, length(simulated.signatures));
          for l=1:nDetected
            Rs4detectedAndSimulatedGeneAxes(l,:) = detectionPCA.PCs(l).corr2simulatedGeneSpacePC;
            Rs4detectedAndSimulatedSampleAxes(l,:) = detectionPCA.PCs(l).corr2simulatedSampleSpacePC;
          end
      end
      comparisons.(sComparisonTarget).Rs4detectedAndSimulatedSampleAxes = Rs4detectedAndSimulatedSampleAxes;
        assert(~any(any(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedSampleAxes>1)), 'some comparisons.(sComparisonTarget).Rs4detectedAndSimulatedSampleAxes are > 1');
      comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes = Rs4detectedAndSimulatedGeneAxes;
        assert(~any(any(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes>1)), 'some comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes are > 1');
    %Extract and save best matches for plotting:
      comparisons.(sComparisonTarget).SI4Presentation = []; %indices relative to .Rs4detectedAndSimulatedGeneAxes and .Rs4detectedAndSimulatedSampleAxes.
      %R4BestmatchSorting = abs(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes) + abs(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedSampleAxes); %average correlaitons in gene and sample space.
      R4BestmatchSorting = abs(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes); %we are plotting correlations in gene space below only.
      [r4MaxI,maxI] = max(R4BestmatchSorting,[],2); %which simulated signature (column maxI) does a gene space PC fit best?
      %Assign each detected signature to that simulated signature to which it has the highest absolute correlation:
        for coli=1:max(maxI) %rows=detected signatures, cols=simulated signatures; add for each col all detected signatures with maximal |r| to this col:
          %which detected signatures (rowI) have highest correlation to the current simulated signature (coli)?
            rowIWithCurrentColAsBestMatch = find(maxI==coli); 
            if(isempty(rowIWithCurrentColAsBestMatch))
              warning('comparison target=%s: no detected signature axis has maximal correlation to simulated signature coli=%d (false negative?)', sComparisonTarget, coli);
            end
            %If there are multiple detected signatures that have maximal correlation to the current simulated signature, sort them descending by their correlation:
              rowIWithCurrentColAsBestMatch_sortedByCorr = ilsub(rowIWithCurrentColAsBestMatch, ilnth(2,@sort,r4MaxI(rowIWithCurrentColAsBestMatch),'descend'));
          %Append current best-matching detected axes to the presentation order of detected signaturers:
            comparisons.(sComparisonTarget).SI4Presentation = [comparisons.(sComparisonTarget).SI4Presentation; rowIWithCurrentColAsBestMatch_sortedByCorr];
        end
        if(length(unique(comparisons.(sComparisonTarget).SI4Presentation))<length(comparisons.(sComparisonTarget).SI4Presentation))
          warning('comparison target=%s: a single detected axis is the best match for more than one simulated signature (this indicates mixing instead of dissection of signatures)', sComparisonTarget);
        end
  end

%% Visualize simulation/detection comparisons:
  fprintf('\nNote: this self test for randomly generated input data can be considered passed, if generated plots resemble corresponding expected results at \\SDCM_selfTest_expectedResults\\selfTest_highDim_exemplaryResults*.png\n');
  %Plot the remaining signal after SDCM dissection:
    if(true)
      f=figure('Position',[73 1 505 997],'Color','w');
      CM_supportsNaNs = [[.5 .5 .5];CM];
      cLim=[-2,2];
      data2Plot = detection.remainingSignal.L2Rs;
        if(true) %support NaN color.
          lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(CM,1);
          data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
        end
      imagesc(data2Plot); 
        caxis(cLim); 
      title(sprintf('leftover noise signal'),'FontSize',8);
      set(gca,'XTick',[],'YTick',[]);
      xlabel(sprintf('%d samples',simulated.nP));
      ylabel(sprintf('%d genes',simulated.nG));
      colormap(CM_supportsNaNs); 
        cbh=colorbar;
        if(true) %remove the gray color from the color bar:
          hImg = findobj(cbh,'Type','Image');
          if(~isempty(hImg)) %colormap image no longer accessible since Matlab's new graphics engine.
            CData = get(hImg,'CData');
            if(CData(end)==size(CM_supportsNaNs,1))
              CData(1)=2;
            end
            set(findobj(cbh,'Type','Image'),'CData',CData);
          end
        end
        daspect([1 1 1]);
    end
  %Plot correlation matrix:
    for sComparisonTarget = {'methodValidation','PCA'}; sComparisonTarget=sComparisonTarget{1}; 
      for sComparisonMetric = {'Rs4detectedAndSimulatedGeneAxes'}; sComparisonMetric=sComparisonMetric{1}; %'Rs4detectedAndSimulatedSampleAxes'
        %Get correlations in presentation order:
          SI = comparisons.(sComparisonTarget).SI4Presentation;
          switch(sComparisonMetric)
            case 'Rs4detectedAndSimulatedSampleAxes'
              M = comparisons.(sComparisonTarget).(sComparisonMetric)(SI,:);
            case 'Rs4detectedAndSimulatedGeneAxes'
              M = comparisons.(sComparisonTarget).(sComparisonMetric)(SI,:);
            case 'mean(Rs4detectedAndSimulatedGeneAxes,Rs4detectedAndSimulatedSampleAxes)'
              M = (...
                 abs(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedGeneAxes(SI,:))...
               + abs(comparisons.(sComparisonTarget).Rs4detectedAndSimulatedSampleAxes(SI,:))...
              )/2;
              M(BIsNotDetected,:) = -Inf;
          end
          %Get the detection order (iteration indices k) respectively the PC rank by variance for each row:
            kRespectivelyPCRanks = cellfun(@num2str,num2cell(SI),'UniformOutput',false); %write detection positions left of the pixel rows.
        %Plot it:
          f=figure('Position',[1+iif(~strcmp(sComparisonTarget,'methodValidation'),685,0), 1, 679, 577],'Color','w');
            data2plot = abs(M);
            imagesc(data2plot); 
            set(gca,'LineWidth',2);
            colormap(CM4Rs);
            caxis([0,1]);
            cbh=colorbar;
            daspect([1 1 1]);
          %Descriptions:
            xlabel('simulated signatures'); 
            ylabel(ils(sComparisonTarget...
              ,'methodValidation',sprintf('detected signatures\n(indices indicate detection iterations k)')...
              ,'PCA',sprintf('detected signatures\n(indices indicate PC ranks by their variance)')...
            ));
            set(gca,'XTick',[],'YTick',1:size(M,1),'YTickLabel',kRespectivelyPCRanks);
            switch(sComparisonMetric)
              case 'Rs4detectedAndSimulatedSampleAxes'
                title(sprintf(...
                  '%s, sample axes correlations%s'...
                 ,sTestMnemonic ...
                 ,ils(sComparisonTarget,'methodValidation',', \bfSDCM results', 'PCA',', \bfPCA results', 'ICA',', \bfICA results', 'NNMF',', \bfNNMF results') ...
                ));
                ylabel(cbh,'sample axes correlations');
              case 'Rs4detectedAndSimulatedGeneAxes'
                title(sprintf(...
                  '%s, gene axes correlations%s'...
                 ,sTestMnemonic ...
                 ,ils(sComparisonTarget,'methodValidation',', \bfSDCM results', 'PCA',', \bfPCA results', 'ICA',', \bfICA results', 'NNMF',', \bfNNMF results') ...
                ));
                ylabel(cbh,'gene axes correlations');
              case 'mean(Rs4detectedAndSimulatedGeneAxes,Rs4detectedAndSimulatedSampleAxes)'
                title(sprintf(...
                  'overview, avg gene and sample axes correlations%s'...
                 ,ils(sComparisonTarget,'methodValidation',', \bfSDCM results', 'PCA',', \bfPCA results', 'ICA',', \bfICA results', 'NNMF',', \bfNNMF results') ...
                ));
                ylabel(sprintf('mean of gene and\nsample axes correlations'));
            end
      end 
    end
