% SYNTAX
%   function [...
%     ISource, ITarget, ...
%     ISourceNonexistentInTarget, ITargetNonexistentInSource, ...
%     ISourceUnexplainedDuplicates, ITargetUnexplainedDuplicates, ...
%     CITarget2Source, CISource2Target ... %enthält den kompletten Kreuzmultijoin
%   ]=IFromTo(...
%     sourceStrings, targetStrings, sPostprocessingModeForDuplicates, bAllowEmptyMatch ...
%    ,sourceAliasesFcn, targetAliasesFcn ... 
%    ,customCompareFcn, customSortCols ...
%   )
% NOTE:
%   The duplicates mode determines the association of source and target indices in 
%   case of ambiguities. Supported duplicate modes are: 'rawJoinArraysOnly', 
%   'matchInOriginalOrderSourceSorting', 'matchInOriginalOrderTargetSorting',
%   'allWithAllSourceSorting', 'allWithAllTargetSorting'.
%   Default is 'matchInOriginalOrderTargetSorting', i.e. if an ID has multiple
%   occurrences in both sourceStrings and targetStrings, the first occurrence
%   in targetStrings is associated to the first occurence in sourceStrings.
%   Henceforth this sorce string will not be associated any more with any other
%   target strings. (This is useful to resolve ID ambiguities via a defined
%   order. If all ambiguities are of same rank (same occurrences of the ambigous
%   ID string in source and target), a one-to-one association is returned nonetheless
%   using the order of IDs to disambiguou.)
% AUTHOR
%  (C) Michael Grau, March 2010.

function [...
      ISource, ITarget ...
     ,ISourceNonexistentInTarget, ITargetNonexistentInSource ...
     ,ISourceUnexplainedDuplicates, ITargetUnexplainedDuplicates ...
     ,CITarget2Source, CISource2Target ...
] = IFromTo(...
      sourceRows, targetRows...
     ,sPostprocessingModeForDuplicates, bAllowEmptyMatch ...
     ,sourceAliasesFcn, targetAliasesFcn ... 
     ,customCompareFcn ...
     ,customSortCols ...
)
  %% Parameter logic and initialization:
    if(nargin<3) 
      sPostprocessingModeForDuplicates = 'matchInOriginalOrderTargetSorting'; 
    else
      knownDuplicateModes = {'rawJoinArraysOnly', 'matchInOriginalOrderSourceSorting', 'matchInOriginalOrderTargetSorting', 'allWithAllSourceSorting', 'allWithAllTargetSorting'};
      assert(ismember(sPostprocessingModeForDuplicates, knownDuplicateModes), 'sPostprocessingModeForDuplicates must be one of the following: %s', cellstring2separatedList(knownDuplicateModes));
    end
    if(nargin<4) 
      bAllowEmptyMatch = false; 
    else
      assert(isscalar(bAllowEmptyMatch) && islogical(bAllowEmptyMatch), 'bAllowEmptyMatch must be a logical scalar');
    end;
    if(nargin<5 || isempty(sourceAliasesFcn)) 
      sourceAliasesFcn = []; 
    else
      assert(isa(targetAliasesFcn,'function_handle'), 'if provided, targetAliasesFcn must be a function handle of signature targetAliases = sourceAliasesFcn(targetSCA)'); 
    end;
    if(nargin<6 || isempty(targetAliasesFcn)) 
      targetAliasesFcn = []; 
    else
      assert(isa(targetAliasesFcn,'function_handle'), 'if provided, targetAliasesFcn must be a function handle of signature targetAliases = sourceAliasesFcn(targetSCA)'); 
    end
    if(nargin<7 || isempty(customCompareFcn)) 
      customCompareFcn = []; 
    else
      assert(isa(customCompareFcn,'function_handle') && (nargin(customCompareFcn)==2 || nargin(customCompareFcn)==3), 'if provided, customCompareFcn must be a function handle of signature @(sourceSCA, targetSCA) or @(sourceSCA, targetSCA, sDomain) returning [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch]'); 
    end
    if(nargin<8 || isempty(customSortCols)) 
      customSortCols = ':';
    else
      assert(strcmp(customSortCols,':') || all(ismember(customSortCols,1:size(sourceRows,2))), 'if provided, customSortCols must be a indes vector with values in 1:size(sourceRows,2)'); 
    end
    bSimpleComparison = isempty(sourceAliasesFcn) && isempty(targetAliasesFcn) && ~isa(customCompareFcn,'function_handle');
    if(size(sourceRows,1)==1 && size(sourceRows,2)>1 || size(targetRows,1)==1 && size(targetRows,2)>1)
      error('parameters sourceRows and targetRows must be provided as column vectors, if there is only one compar dimension; else they will be treated like one row with length(sourceRows) compare dimensions.');
    end
  %% define and select comparer functions:
    %simple comparer for numeric types:
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_forNumericTypes_singleComparisonDim(sourceSCA, targetSCA)
        n1 = sourceSCA{1};
        n2 = targetSCA{1};
        if(isempty(n1) && isempty(n2))
          bIsEmptyMatch = true;
          bIsEqual = true;
          bSourceSCA_gte_TargetSCA = true;
        else
          bIsEmptyMatch = false;
          bIsEqual = n1==n2;
          bSourceSCA_gte_TargetSCA = n1>=n2;
        end
      end    
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_forNumericTypes_singleComparisonDim_scalar(sourceSCA, targetSCA)
        bIsEmptyMatch = false;
        bIsEqual = sourceSCA==targetSCA;
        bSourceSCA_gte_TargetSCA = sourceSCA>=targetSCA;
      end    
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_forNumericTypes_multiComparisonDims(sourceSCA, targetSCA)
        %fast comparison method, if there are no alias functions to use:
          bIsEqual = true;
          bSourceSCA_gte_TargetSCA = true;
          bIsEmptyMatch = true;
          if(iscell(sourceSCA) && iscell(targetSCA))
            for subj=1:length(sourceSCA)
              if(isempty(sourceSCA{subj}))
                if(isempty(targetSCA{subj}))
                  %b still true
                  %bSourceSCA_gte_TargetSCA still true
                  %bIsEmptyMatch still true;
                  continue; %go to next comparison dimension of source and target cell rows.
                end
              end
              bIsEmptyMatch = false;
              bIsEqual = bIsEqual && sourceSCA{subj} == targetSCA{subj};
              if(~bIsEqual) %bSourceSCA_gte_TargetSCA is only needed, if ~bIsEqual.
                bSourceSCA_gte_TargetSCA = sourceSCA{subj} >= targetSCA{subj};
                break;
              end
            end
          elseif(isnumeric(sourceSCA) && isnumeric(targetSCA))
            bIsEmptyMatch = false;
            for subj=1:length(sourceSCA)
              bIsEqual = bIsEqual && sourceSCA(subj) == targetSCA(subj);
              if(~bIsEqual) %bSourceSCA_gte_TargetSCA is only needed, if ~bIsEqual.
                bSourceSCA_gte_TargetSCA = sourceSCA(subj) >= targetSCA(subj);
                break;
              end
            end
          end
      end
    %simple comparer for string types:
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_forStringTypes_multiComparisonDims(sourceSCA, targetSCA)
        %fast comparison method, if there are no alias functions to use:
          bIsEqual = true;
          bSourceSCA_gte_TargetSCA = NaN; %true;
          bIsEmptyMatch = true;
          for subj=1:length(sourceSCA)
            if(isempty(sourceSCA{subj}))
              if(isempty(targetSCA{subj}))
                %b still true
                %bSourceSCA_gte_TargetSCA still true
                %bIsEmptyMatch still true;
                continue; %go to next comparison dimension of source and target cell rows.
              end
            end
            bIsEmptyMatch = false;
            bIsEqual = bIsEqual && strcmp(sourceSCA{subj}, targetSCA{subj});
            if(~bIsEqual) %bSourceSCA_gte_TargetSCA is only needed, if ~bIsEqual.
              if(nargout>=2) %Performance: only compute relative string order, if needed:
                minLength = min(length(sourceSCA{subj}), length(targetSCA{subj}));
                letterComparisonResults = sign(sourceSCA{subj}(1:minLength)-targetSCA{subj}(1:minLength));
                signOfFirstDifferentLetter = letterComparisonResults(find(letterComparisonResults,1,'first'));
                if(isempty(signOfFirstDifferentLetter)) %all minLength letters are equal; in this case the shorter string is smaller ('b'<'bb' and ''<'a'):
                  bSourceSCA_gte_TargetSCA = length(sourceSCA{subj}) > length(targetSCA{subj});
                else
                  bSourceSCA_gte_TargetSCA = signOfFirstDifferentLetter==1;
                end
              end
              break;
            end
          end
      end    
    %simple comparer for mixed numeric/string cell arrays:
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_forMixedNumericStringArrays_multiComparisonDims(sourceSCA, targetSCA)
        %fast comparison method, if there are no alias functions to use:
          bIsEqual = true;
          bSourceSCA_gte_TargetSCA = true;
          bIsEmptyMatch = true;
          for subj=1:length(sourceSCA)
            if(isempty(sourceSCA{subj}))
              if(isempty(targetSCA{subj}))
                %b still true
                %bSourceSCA_gte_TargetSCA still true
                %bIsEmptyMatch still true;
                continue; %go to next comparison dimension of source and target cell rows.
              end
            end
            bIsEmptyMatch = false;
            if(isnumeric(sourceSCA{subj}) && isnumeric(targetSCA{subj}))
              bIsEqual = bIsEqual && sourceSCA{subj} == targetSCA{subj};
              if(~bIsEqual) %bSourceSCA_gte_TargetSCA is only needed, if ~bIsEqual.
                bSourceSCA_gte_TargetSCA = sourceSCA{subj} >= targetSCA{subj};
                break;
              end
            elseif(ischar(sourceSCA{subj}) && ischar(targetSCA{subj}))
              bIsEqual = bIsEqual && strcmp(sourceSCA{subj}, targetSCA{subj});
              if(~bIsEqual) %bSourceSCA_gte_TargetSCA is only needed, if ~bIsEqual.
                if(nargout>=2) %Performance: only compute relative string order, if needed:
                  minLength = min(length(sourceSCA{subj}), length(targetSCA{subj}));
                  letterComparisonResults = sign(sourceSCA{subj}(1:minLength)-targetSCA{subj}(1:minLength));
                  signOfFirstDifferentLetter = letterComparisonResults(find(letterComparisonResults,1,'first'));
                  if(isempty(signOfFirstDifferentLetter)) %all minLength letters are equal; in this case the shorter string is smaller ('b'<'bb' and ''<'a'):
                    bSourceSCA_gte_TargetSCA = length(sourceSCA{subj}) > length(targetSCA{subj});
                  else
                    bSourceSCA_gte_TargetSCA = signOfFirstDifferentLetter==1;
                  end
                end
                break;
              end
            else
              error('cannot compare numbers with strings');
            end
          end
      end
    %comparer calculating equality based on sourceAliasesFcn and targetAliasesFcn:
      function [bIsEqual, bSourceSCA_gte_TargetSCA, bIsEmptyMatch] = checkEquality_usingSourceAndTargetAliasFunctions(sourceSCA, targetSCA)
error('checkEquality_usingSourceAndTargetAliasFunctions needs a recheck on example/unit test; especially wrt. bSourceSCA_gte_TargetSCA handling in this case');
        %parameter check:
          if(size(sourceSCA,1)~=1 || size(targetSCA,1)~=1) error('internal: comparison must be row-wise'); end;
          if(size(sourceSCA,2)~=size(targetSCA,2)) error('internal: inconsistent number of source and target string fields'); end;
          bIsEmptyMatch = true;
        %apply alias functions:
          sourceAliases = cellfun(@(S){S},sourceSCA,'UniformOutput',false);
            if(~isempty(sourceAliasesFcn)) 
              sourceAliases = sourceAliasesFcn(sourceSCA); %must return all aliases to sourceSCA(1,i) in sourceAliases{i} as cellstring array, including sourceSCA(1,i).
              bIsEmptyMatch = false;
              if(length(sourceAliases)~=size(sourceSCA,2)) error('internal: alias function must return a cell of cells, where aliases{j} contains a cellstring array for all aliases to sourceSCA{j}, under consideration of all columns sourceSCA{:}'); end;
            end; 
          targetAliases = cellfun(@(S){S},targetSCA,'UniformOutput',false);
            if(~isempty(targetAliasesFcn)) 
              targetAliases = targetAliasesFcn(targetSCA); %must return all aliases to targetSCA(1,i) in targetAliases{i} as cellstring array, including targetSCA(1,i).
              bIsEmptyMatch = false;
              if(length(targetAliases)~=size(sourceSCA,2)) error('internal: alias function must return a cell of cells, where aliases{j} contains a cellstring array for all aliases to targetSCA{j}, under consideration of all columns targetSCA{:}'); end;
            end;    
        %in every column there must be at least one matching alias:
          bIsEqual = true;
          bSourceSCA_gte_TargetSCA = true;
          for subj=1:size(sourceSCA,2)
            bIsEqual = bIsEqual && ~isempty(intersect(sourceAliases{subj}, targetAliases{subj}));
            bSourceSCA_gte_TargetSCA = bSourceSCA_gte_TargetSCA && length(sourceAliases{subj})>=length(targetAliases{subj});
          end
      end    
    %assign comparer function for this call:
      if(bSimpleComparison)
        if(  (isnumeric(sourceRows) || iscell(sourceRows) && all(cellfun(@isnumeric,sourceRows(:))))...
          && (isnumeric(targetRows) || iscell(targetRows) && all(cellfun(@isnumeric,targetRows(:))))...
        )
          if(size(sourceRows,2)==1 && size(targetRows,2)==1) %only one comparison dimension:
            %Performance: if we have only numeric arrays, skip empty checks:
              bScalarNonEmptyIDs = isnumeric(sourceRows) && isnumeric(targetRows) || all(cellfun(@isscalar,sourceRows)) && all(cellfun(@isscalar,targetRows));
              if(bScalarNonEmptyIDs)
                if(~isnumeric(sourceRows))
                  sourceRows = vertcat(sourceRows{:});
                end
                if(~isnumeric(targetRows))
                  targetRows = vertcat(targetRows{:});
                end
                checkEquality = @checkEquality_forNumericTypes_singleComparisonDim_scalar;
            %Standard case supporting empty cells:
              else
                checkEquality = @checkEquality_forNumericTypes_singleComparisonDim;
              end
          else %>1 comparison dimensions:
            checkEquality = @checkEquality_forNumericTypes_multiComparisonDims;
          end
        elseif(  iscellstr(sourceRows) && iscellstr(targetRows))
          checkEquality = @checkEquality_forStringTypes_multiComparisonDims;
        elseif(all(cellfun(@(c)isnumeric(c)||ischar(c), sourceRows(:))) && all(cellfun(@(c)isnumeric(c)||ischar(c), targetRows(:))))
          checkEquality = @checkEquality_forMixedNumericStringArrays_multiComparisonDims;
        else
          error('input data were not basic (either all numbers or all strings), but neither customCompareFcn nor sourceAliasesFcn/targetAliasesFcn were specified.');
        end
      elseif(isa(customCompareFcn,'function_handle'))
        checkEquality = customCompareFcn;
      elseif(~isempty(sourceAliasesFcn) || ~isempty(targetAliasesFcn))
        checkEquality = @checkEquality_usingSourceAndTargetAliasFunctions;
      else
        error('unknown comparer function; see help.');
      end
      nArgin4checkEquality = nargin(checkEquality);

  %% Temporary internal sorting for comparison:
    bUseAllCols4Sorting = strcmp(customSortCols, ':');
    if(bUseAllCols4Sorting)
      [SsourceRows, SI1] = sortrows(sourceRows, 1:size(sourceRows,2));
      [StargetRows, SI2] = sortrows(targetRows, 1:size(targetRows,2));
    else
      [SsourceRows, SI1] = sortrows(sourceRows, iif(ischar(customSortCols),@()ilsub(1:size(sourceRows,2),customSortCols),customSortCols));
      [StargetRows, SI2] = sortrows(targetRows, iif(ischar(customSortCols),@()ilsub(1:size(targetRows,2),customSortCols),customSortCols));
    end
  %% finde die pointer cursorSource->cursorTarget und speichere in ITarget:
    cursorSource = 1; cursorTarget = 1;
    CISource2Target = cell(size(SsourceRows,1),1); %CISource2Target(i) enth�lt alle target indizes, die SsourceStrings{i,:} gleichen
    CITarget2Source = cell(size(StargetRows,1),1); %CITarget2Source(j) enth�lt alle source indizes, die StargetStrings{j,:} gleichen
    lastTime = clock; lastTime = lastTime(end); lastCursorSum=2;
    while(cursorSource <= size(SsourceRows,1) && cursorTarget <= size(StargetRows,1))
      %status info:
        if(cursorSource+cursorTarget > lastCursorSum+ceil(size(SsourceRows,1)/1000))
          lastCursorSum = cursorSource+cursorTarget;
          currentSec = clock; currentSec = currentSec(end); 
          if(mod(-lastTime+currentSec, 60) > 25)
            prc = 100*max(cursorSource/size(SsourceRows,1), cursorTarget/size(StargetRows,1));
            disp(['Index calculation IFromTo / ',num2str(prc,'%0.1f'),'%']); drawnow;
            lastTime = clock; lastTime = lastTime(end); 
          end
        end
      %calculate equality and sort order of SsourceStrings(cursorSource,:), StargetStrings(cursorTarget,:):
        if(nArgin4checkEquality==2)
          [bIsEqual, bSourceRow_gte_TargetRow, bIsEmptyMatch] = checkEquality(SsourceRows(cursorSource,:), StargetRows(cursorTarget,:));
        else %support for comparers that need to know the compared domains:
          [bIsEqual, bSourceRow_gte_TargetRow, bIsEmptyMatch] = checkEquality(SsourceRows(cursorSource,:), StargetRows(cursorTarget,:), 'sourceVsTarget');
        end
        %bEmptyMatchAllowedOrNotEmpty = bAllowEmptyMatch || (~isempty([SsourceStrings{cursorSource,:}]) && (~all(isnumeric([SsourceStrings{cursorSource,:}])) || ~any(isnan([SsourceStrings{cursorSource,:}]))));
      %falls match:
        if(bIsEqual && (bAllowEmptyMatch || ~bIsEmptyMatch))
          %target duplicates lookahead:
            targetLookahead = 0;
            %support range comparers/Apr2014: if a target range element matches the same source elements, we must not advance the target range element until the next source element is strictly greater than the current target range element:
              bMayAdvanceTarget = true;
              if(~isempty(customCompareFcn) && cursorSource+1<=size(SsourceRows,1)) %interne comparer sind keine range comparer && wenn die Quelle noch nicht ans Ende gespult ist, teste:
                if(nArgin4checkEquality==2)
                  [bIsEqualWithNextSourceRow, bNextSourceRow_gte_TargetRow, bIsEmptyMatch2] = checkEquality(SsourceRows(cursorSource+1,:), StargetRows(cursorTarget,:));
                else %support for comparers that need to know the compared domains:
                  [bIsEqualWithNextSourceRow, bNextSourceRow_gte_TargetRow, bIsEmptyMatch2] = checkEquality(SsourceRows(cursorSource+1,:), StargetRows(cursorTarget,:), 'sourceVsTarget');
                end
                bMayAdvanceTarget = bMayAdvanceTarget && ~bIsEqualWithNextSourceRow; %die n�chste Quellzeile darf nicht in der Range der aktuellen Zielzeile liegen
                bMayAdvanceTarget = bMayAdvanceTarget && bNextSourceRow_gte_TargetRow; %die n�chste Quellzeile muss strikt gr��er als die aktuelle Zielzeile sein
              end
            if(bMayAdvanceTarget)
              targetLookahead = 1;
              while(true)
                if(nArgin4checkEquality==2)
                  bContinue = cursorTarget+targetLookahead-1 < size(StargetRows,1) && checkEquality(StargetRows(cursorTarget+targetLookahead,:),StargetRows(cursorTarget,:));
                else %support for comparers that need to know the compared domains:
                  bContinue = cursorTarget+targetLookahead-1 < size(StargetRows,1) && checkEquality(StargetRows(cursorTarget+targetLookahead,:),StargetRows(cursorTarget,:),'targetVsTarget');
                end
                if(~bContinue) break; end
                %<-KORREKTUR am 10.09.2012: "-1 < ..."; vorher wurde das letzte Element, falls es ein Duplikat war, nicht korrekt zugeordnet!
                targetLookahead = targetLookahead + 1;
              end
            end
            ITargetDuplicates = cursorTarget : cursorTarget+max(1,targetLookahead)-1; %in case of target range elements, targetLookahead may be zero, but we want to associate the range element to the current source element, so max(1,.).
          %source duplicates lookahead:
            sourceLookahead = 0;
            %support range comparers/Apr2014: if a source range element matches the same target elements, we must not advance the source range element until the next target element is strictly greater than the current source range element:
              bMayAdvanceSource = true;
              if(~isempty(customCompareFcn) && cursorTarget+1<=size(StargetRows,1)) %interne comparer sind keine range comparer && wenn das Ziel noch nicht ans Ende gespult ist, teste:
                if(nArgin4checkEquality==2)
                  [bIsEqualWithNextTargetRow, bSourceRow_gte_NextTargetRow, bIsEmptyMatch2] = checkEquality(SsourceRows(cursorSource,:), StargetRows(cursorTarget+1,:));
                else %support for comparers that need to know the compared domains:
                  [bIsEqualWithNextTargetRow, bSourceRow_gte_NextTargetRow, bIsEmptyMatch2] = checkEquality(SsourceRows(cursorSource,:), StargetRows(cursorTarget+1,:), 'sourceVsTarget');
                end
                bMayAdvanceSource = bMayAdvanceSource && ~bIsEqualWithNextTargetRow; %die n�chste Zielzeile darf nicht in der Range der aktuellen Quellzeile liegen
                bMayAdvanceSource = bMayAdvanceSource && ~bSourceRow_gte_NextTargetRow; %die n�chste Zielzeile muss strikt gr��er sein als die aktuelle Quellzeile
              end
            if(bMayAdvanceSource)
              sourceLookahead = 1;
              while(true)
                if(nArgin4checkEquality==2)
                  bContinue = cursorSource+sourceLookahead-1 < size(SsourceRows,1) && checkEquality(SsourceRows(cursorSource+sourceLookahead,:),SsourceRows(cursorSource,:));
                else %support for comparers that need to know the compared domains:
                  bContinue = cursorSource+sourceLookahead-1 < size(SsourceRows,1) && checkEquality(SsourceRows(cursorSource+sourceLookahead,:),SsourceRows(cursorSource,:), 'sourceVsSource');
                end
                if(~bContinue) break; end
                %<-KORREKTUR am 10.09.2012: "-1 < ..."; vorher wurde das letzte Element, falls es ein Duplikat war, nicht korrekt zugeordnet!
                sourceLookahead = sourceLookahead + 1;
              end
            end
            ISourceDuplicates = cursorSource : cursorSource+max(1,sourceLookahead)-1; %in case of source range elements, sourceLookahead may be zero, but we want to associate the source element to the current target element, so max(1,.).
          %speichere match(es):
            %[CISource2Target{ISourceDuplicates}] = deal(ITargetDuplicates);
            %[CITarget2Source{ITargetDuplicates}] = deal(ISourceDuplicates);
            CISource2Target(ISourceDuplicates) = {ITargetDuplicates};
            CITarget2Source(ITargetDuplicates) = {ISourceDuplicates};
          %Advance on the sorted bands:
            cursorTarget = cursorTarget + targetLookahead;
            cursorSource = cursorSource + sourceLookahead;
      %falls ungleich:
        else
          %berechne relative lexikografische Position:
            %[~, relI] = sortrows([SsourceStrings(cursorSource,:); StargetStrings(cursorTarget,:)]);
          %2017-06 hotfix: if ~bSourceRow_gte_TargetRow and we have bIsEmptyMatch, we can advance both the source row or the target row:
            if(~bAllowEmptyMatch && bIsEmptyMatch)
              cursorSource = cursorSource + 1;
              cursorTarget = cursorTarget + 1;
              continue;
            end
          %falls source<target in Sortierreihenfolge, existiert source nicht in target:
            %if(relI(1)==1)
            if(~bSourceRow_gte_TargetRow)
              %CISource2Target(cursorSource) = [];
              cursorSource = cursorSource + 1;
          %falls source>target in Sortierreihenfolge, existiert target nicht in source:
            else
              %CITarget2Source(cursorTarget) = [];
              cursorTarget = cursorTarget + 1;
            end
        end
    end
    %cursorSource
    %cursorTarget
  %% Undo temporary internal sorting:
    %In den Zellen von CITarget2Source stehen Indizes zu SsourceStrings. Wir wollen die korrespondierenden Indizes im unsortierten Array sourceStrings; diese stehen gerede in SI1:
      CITarget2Source = cellfun(@(I)SI1(I), CITarget2Source, 'UniformOutput', false);
      bCISource2TargetNeeded = nargout>=8 || ~strcmp(sPostprocessingModeForDuplicates,'rawJoinArraysOnly');
      if(bCISource2TargetNeeded) %Performance: the following line was profiled to be slow, but is never needed for a call by joinAndFill.
        CISource2Target = cellfun(@(I)SI2(I), CISource2Target, 'UniformOutput', false);
      else
        clear CISource2Target;
      end
    %Die Zellen 1:length(CITarget2Source) von CITarget2Source selbst korrespondieren mit den Eintr�gen von StargetStrings. Wir wollen CITarget2Source so umsortieren, dass danach die Zellen mit den Eintr�gen des originalen Arrays targetStrings korrespondieren. Also z.B. soll die erste Zelle CITarget2Source{1} danach an der Position stehen, wo StargetStrings{1} im originalen Array targetStrings stand. Da targetStrings{SI2}==StargetStrings, stand StargetStrings{1} an SI2(1)-ter Stelle. Es muss also gelten: CISourceNeu{SI2(1)}:=CISourceAlt{1} oder allgemeiner: CISourceNeu{SI2}:=CISourceAlt:
      CITarget2Source(SI2) = CITarget2Source;
      if(bCISource2TargetNeeded)
        CISource2Target(SI1) = CISource2Target;
      end

  %% Postprocessing/erstelle Verbindung und behandele Duplikate:
    %Early break, if the caller only needs the join and no postprocessing for duplicates:
      if(strcmp(sPostprocessingModeForDuplicates,'rawJoinArraysOnly'))
        ISource=[]; ITarget=[];
        ISourceNonexistentInTarget=[]; ITargetNonexistentInSource=[];
        ISourceUnexplainedDuplicates=[]; ITargetUnexplainedDuplicates=[];
        return;
      end
    %Extrahiere unerklärte Elemente:
      ITargetNonexistentInSource = find(cellfun(@isempty, CITarget2Source));
      ISourceNonexistentInTarget = find(cellfun(@isempty, CISource2Target));
    %Logic for sPostprocessingModeForDuplicates 'one2one':
      if(strcmp(sPostprocessingModeForDuplicates,'one2one'))
        sPostprocessingModeForDuplicates = 'matchInOriginalOrderTargetSorting';
        bErrorOnDuplicatesOrMissing = true;
      else
        bErrorOnDuplicatesOrMissing = false;
      end
    %Process duplicates:
      switch(sPostprocessingModeForDuplicates)
        case 'matchInOriginalOrderSourceSorting' 
          %inkrementelles Entfernen erkl�rter targets:
            CITargetCopy = CISource2Target;
            IexplainedTargets = [];
            for i=1:size(sourceRows,1)
              CITargetCopy{i} = setdiff(CITargetCopy{i}, IexplainedTargets);
              IexplainedTargets = [IexplainedTargets; CITargetCopy{i}(1:min(1,end))]; %#ok<AGROW>
            end
          ISource = find(~cellfun(@isempty, CITargetCopy)); %zu welchem source entry existiert mind. ein target?
          ITarget = cellfun(@(I) I(1:min(1,length(I))), CITargetCopy, 'UniformOutput', false); %nehme immer das erste noch nicht erkl�rte target.
            ITarget = vertcat(ITarget{:});
          ISourceUnexplainedDuplicates = setdiff(find(cellfun(@isempty, CITargetCopy)), ISourceNonexistentInTarget);
          ITargetUnexplainedDuplicates = cellfun(@(I) I(2:length(I)), CITargetCopy, 'UniformOutput', false);
            ITargetUnexplainedDuplicates = setdiff(unique(vertcat(ITargetUnexplainedDuplicates{:})), ITarget);          
          %assert(all(sourceRows(ISource)==targetRows(ITarget)))
        case 'matchInOriginalOrderTargetSorting' 
          %inkrementelles Entfernen bereits erkl�rter targets (f�r reihenfolgendefinierte ID-Zusammengeh�rigtkeit bei Anwesenheit von Duplikaten):
            CISourceCopy = CITarget2Source;
            IexplainedSources = zeros(size(sourceRows,1),1); iCursor = 1;
            for i=1:size(targetRows,1)
              if(isempty(CISourceCopy{i}))
                %nothing to do.
              elseif(length(CISourceCopy{i})==1) %exactly one source ID matches the ith target.
                if(ismember(CISourceCopy{i}, IexplainedSources))
                  CISourceCopy{i} = [];
                end
              else %multiple source IDs match the ith target; remove those that already have previos targets matching it in original order (target sorting).
                CISourceCopy{i} = setdiff(CISourceCopy{i}, IexplainedSources);
              end
              if(~isempty(CISourceCopy{i}))
                %IexplainedSources = [IexplainedSources, CISourceCopy{i}(1)];
                IexplainedSources(iCursor) = CISourceCopy{i}(1);
                iCursor = iCursor + 1;
              end
              if(mod(i,10000)==0)
                disp([num2str(100*i/size(targetRows,1)),'% postprocessed']);
              end
            end
          ITarget = find(~cellfun(@isempty, CISourceCopy));
          ISource = cellfun(@(I) I(1:min(1,length(I))), CISourceCopy, 'UniformOutput', false);
            ISource = vertcat(ISource{:});
          ITargetUnexplainedDuplicates = setdiff(find(cellfun(@isempty, CISourceCopy)), ITargetNonexistentInSource);
          ISourceUnexplainedDuplicates = cellfun(@(I) reshape(I(2:length(I)),1,max(0,length(I)-1)), CISourceCopy, 'UniformOutput', false);
            ISourceUnexplainedDuplicates = setdiff(unique([ISourceUnexplainedDuplicates{:}]), ISource);
          %assert(all(sourceRows(ISource)==targetRows(ITarget)))
        case 'allWithAllSourceSorting'
          %CTemp = cellfun(@(I)I,CISource2Target,'UniformOutput',false);
          ITarget = vertcat(CISource2Target{:});
          L = cellfun(@length, CISource2Target);
          ISource = nan(sum(L),1);
          cursor = 1;
          for i=1:size(sourceRows,1)
           ISource(cursor:cursor-1+L(i)) = i; %March2015: added -1
           cursor = cursor+L(i);
          end
          %ISource = [];
          %for i=1:size(sourceStrings,1)
          %  ISource = [ISource, i*ones(1,L(i))];
          %end
          ISourceUnexplainedDuplicates = [];
          ITargetUnexplainedDuplicates = [];
          %assert(all(sourceRows(ISource)==targetRows(ITarget)))
        case 'allWithAllTargetSorting'
          %CTemp = cellfun(@(I)I,CITarget2Source,'UniformOutput',false);
          ISource = vertcat(CITarget2Source{:});
          L = cellfun(@length, CITarget2Source);
          ITarget = nan(sum(L),1);
          cursor = 1;
          for i=1:size(targetRows,1)
           ITarget(cursor:cursor-1+L(i)) = i; %March2015: added -1
           cursor = cursor+L(i);
          end
          %ITarget = [];
          %for i=1:size(targetStrings,1)
          %  ITarget = [ITarget, i*ones(1,L(i))];
          %end        
          ITargetUnexplainedDuplicates = [];
          ISourceUnexplainedDuplicates = [];
          %assert(all(sourceRows(ISource)==targetRows(ITarget)))
        otherwise 
          error('Unknown posprocessing mode for duplicates "%s"; see help.', sPostprocessingModeForDuplicates);
      end
    %Error/validation logic, if enabled:
      if(bErrorOnDuplicatesOrMissing)
        if(~isempty(ISourceNonexistentInTarget))
          error(['not all source strings are available in the target: ', sourceRows{ISourceNonexistentInTarget,:}]);
        end
        if(~isempty(ITargetNonexistentInSource))
          error(['not all target strings are available in the source: ', sourceRows{ITargetNonexistentInSource,:}]);
        end
        if(~isempty(ITargetUnexplainedDuplicates))
          error(['some target strings are not unique and do not have a corresponding source duplicate: ', sourceRows{ITargetUnexplainedDuplicates,:}]);
        end
        if(~isempty(ISourceUnexplainedDuplicates))
          error(['some source strings are not unique and do not have a corresponding target duplicate: ', sourceRows{ISourceUnexplainedDuplicates,:}]);
        end
      end
end

%% Examples: 
%     clc
%     sourceStrings = {'1','1','2','3','2','7','',''}
%     targetStrings = {'','1','3','2','1','1',''}
%     bAllowEmptyMatch = false;
%     sPostprocessingModeForDuplicates = 'allWithAllTargetSorting';
% %     sPostprocessingModeForDuplicates = 'matchInOriginalOrderSourceSorting';    
% %     sPostprocessingModeForDuplicates = 'matchInOriginalOrderTargetSorting';    
%     [ ISource, ITarget, ...
%       ISourceNonexistentInTarget, ITargetNonexistentInSource, ...
%       ISourceUnexplainedDuplicates, ITargetUnexplainedDuplicates, ...
%       CITarget2Source, CISource2Target ... %enth�lt den kompletten Kreuzmultijoin
%     ] = IFromTo(...
%       sourceStrings, targetStrings, bAllowEmptyMatch, sPostprocessingModeForDuplicates ...
%     );
%     %Ausgabe:
%       disp('sourceStrings(ISource): '); disp(sourceStrings(ISource)); %source data in target order
%       disp('targetStrings(ITarget): '); disp(targetStrings(ITarget)); %joined data from target set
%       disp('sourceStrings(ISourceNonexistentInTarget): '); disp(sourceStrings(ISourceNonexistentInTarget));
%       disp('targetStrings(ITargetNonexistentInSource): '); disp(targetStrings(ITargetNonexistentInSource));
%       disp('sourceStrings(ISourceUnexplainedDuplicates): '); disp(sourceStrings(ISourceUnexplainedDuplicates));
%       disp('targetStrings(ITargetUnexplainedDuplicates): '); disp(targetStrings(ITargetUnexplainedDuplicates));

%     clc
%     sourceStrings = {DB.Signatures.NFkappaB.probesMetaInfo.GlobalAcc}';
%     targetStrings = {DB.probesMetaInfo.GlobalAcc}';
%     bAllowEmptyMatch = false;
%     sPostprocessingModeForDuplicates = 'allWithAllTargetSorting';
%     sPostprocessingModeForDuplicates = 'matchInOriginalOrderSourceSorting';    
%     sPostprocessingModeForDuplicates = 'matchInOriginalOrderTargetSorting';    
%     [ ISource, ITarget, ...
%       ISourceNonexistentInTarget, ITargetNonexistentInSource, ...
%       ISourceUnexplainedDuplicates, ITargetUnexplainedDuplicates, ...
%       CITarget2Source, CISource2Target ... %enth�lt den kompletten Kreuzmultijoin
%     ] = IFromTo(...
%       sourceStrings, targetStrings, bAllowEmptyMatch, sPostprocessingModeForDuplicates ...
%     );
%     %Ausgabe:
%       %disp('sourceStrings(ISource): '); disp(sourceStrings(ISource)); %source data in target order
%       %disp('targetStrings(ITarget): '); disp(targetStrings(ITarget)); %joined data from target set
%       disp('sourceStrings(ISourceNonexistentInTarget): '); disp(sourceStrings(ISourceNonexistentInTarget));
%       %disp('targetStrings(ITargetNonexistentInSource): '); disp(targetStrings(ITargetNonexistentInSource));
%       disp('sourceStrings(ISourceUnexplainedDuplicates): ');
%         for i=ISourceUnexplainedDuplicates'
%           disp(['  sourceStrings{',mat2str(CITarget2Source{CISource2Target{i}}),'} all== ',sourceStrings{i}, ', ',num2str(i),' left unexplained.']);
%         end
%       disp('targetStrings(ITargetUnexplainedDuplicates): ');
%         for i=ITargetUnexplainedDuplicates'
%           disp(['  targetStrings{',mat2str(CISource2Target{CITarget2Source{i}}),'} all== ',targetStrings{i}, ', ',num2str(i),' left unexplained.']);
%         end

%% Code backup:
%   %initalize:
%     if(nargin<3) bAllowEmptyMatch = false; end;
%   %tempor�re interne Sortierung:
%     [SsourceStrings, SI1] = sort(sourceStrings);
%     [StargetStrings, SI2] = sort(targetStrings);
%   %finde die pointer cursorSource->cursorTarget und speichere in ITarget:
%     cursorSource = 1; cursorTarget = 1;
%     ITarget = zeros(size(SsourceStrings,1),1);
%     ITargetUnjoinedDuplicates = [];
%     SAmbigiousTargets = {}; SAmbigiousSources = {};
%     while(cursorSource<=size(SsourceStrings,1))
%       %empty strings: keine matches von empty strings durchf�hren (werden unten wie "NotInTarget" gewertet):
%         if(~bAllowEmptyMatch && isempty(SsourceStrings{cursorSource,:}))
%           ITarget(cursorSource) = -1; %=>ISourceNotInTarget.
%           cursorSource = cursorSource + 1;
%       %falls match:
%         elseif(strcmp(SsourceStrings(cursorSource), StargetStrings(cursorTarget)))
%           %speichere match:
%             ITarget(cursorSource) = cursorTarget;
%           %target duplicates: checke, ob mehrfach in "unique"StringSuperset; falls ja, merke in ITargetUnjoinedDuplicates:
%             %aufr�umen: entferne aus ITargetUnjoinedDuplicates, falls jetzt auch in sourceStrings ein Duplikat war, das dieses erkl�ren kann:
%               previousLookaheadTargetDuplicateI = find(ITargetUnjoinedDuplicates==cursorTarget);
%               if(~isempty(previousLookaheadTargetDuplicateI))
%                 SAmbigiousTargets{end+1} = StargetStrings{ITargetUnjoinedDuplicates(previousLookaheadTargetDuplicateI)};
%                 SAmbigiousSources{end+1} = SsourceStrings{cursorSource,:};
%                 ITargetUnjoinedDuplicates(previousLookaheadTargetDuplicateI) = [];
%               end
%             %lookahead �ber alle Duplikate in "unique"StringSuperset:
%               lookahead = 1;
%               while(cursorTarget+lookahead<length(targetStrings) && strcmp(targetStrings{cursorTarget+lookahead},targetStrings{cursorTarget,:}))
%                 lookahead = lookahead + 1;
%               end
%               ITargetLocalDuplicates = cursorTarget+1 : cursorTarget+lookahead-1;
%               ITargetUnjoinedDuplicates = unique([ITargetUnjoinedDuplicates, ITargetLocalDuplicates]);
%             %cursorTarget um eins inkrementieren, falls es Duplikate gibt:
%               if(~isempty(ITargetLocalDuplicates))
%                 cursorTarget = cursorTarget + 1;
%               end
%           cursorSource = cursorSource + 1;
%       %falls ungleich:
%         else
%           %berechne relative lexikografische Position:
%             [dummy, relI] = sort([SsourceStrings(cursorSource), StargetStrings(cursorTarget)]);
%           %falls source<target in Sortierreihenfolge, existiert source nicht in target => markieren und cursorTarget nicht hochfahren:
%             if(relI(1)==1)
%               ITarget(cursorSource) = -1; %=>ISourceNotInTarget.
%               cursorSource = cursorSource + 1;
%           %ansonsten (source>target in Sortierreihenfolge) => cursorTarget hochfahren:
%             else
%               cursorTarget = cursorTarget + 1;
%             end
%         end
%     end    
%   %unmatched extrahieren und entfernen:
%     ISource = 1:length(ITarget);
%     ISourceNotInTarget = find(ITarget==-1);
%     ITarget(ISourceNotInTarget) = [];
%     ISource(ISourceNotInTarget) = [];
%     ITargetNotInSource = setdiff(1:size(StargetStrings,1), ITarget);
%   %Duplikate entfernen:
%     IDuplicates = diff([0;ITarget(:)])==0; %ITarget is always increasing here.
%     ISourceUnjoinedDuplicates = ISource(IDuplicates);
%     ISource(IDuplicates) = [];
%     ITarget(IDuplicates) = [];
%   %interne Sortierung zum schnellen Vergleichen r�ckg�ngig machen:
%     ISource = SI1(ISource);
%     ISourceNotInTarget = SI1(ISourceNotInTarget);
%     ISourceUnjoinedDuplicates = SI1(ISourceUnjoinedDuplicates);
%     ITarget = SI2(ITarget);    
%     ITargetUnjoinedDuplicates = SI2(ITargetUnjoinedDuplicates);
%     ITargetNotInSource = SI2(ITargetNotInSource);
%   %Reihenfolge von targetStrings annehmen:
%     [ITarget, SITarget] = sort(ITarget);
%     ISource = ISource(SITarget);
%example: 
%     clc
%     C1 = {'1','1','3','3','7','5','4'}
%     C2 = fliplr({'5','6','7','8','9','1','2','3','4'})
%     [ISource, ITarget] = ISourceTo(C1, C2)
%     C1(ISource) %duplicate free data in target order
%     C2(ITarget) %available data from target set    

%     clc
%     C1 = {rawdata.HLB1.IKKbeta(1).probesMetaInfo(1:3).AgilentFeatureIDs}'
%     C2 = {DB.probesMetaInfo.AgilentFeatureIDs}';
%     [ISource, ITarget] = ISourceTo(C1, C2)
%     C1(ISource) %duplicate free data in target order
%     C2(ITarget) %available data from target set    

%     clc
%     C1 = {DB.Signatures.NFkappaB.probesMetaInfo.GlobalAcc}';
%     C2 = {DB.probesMetaInfo.GlobalAcc}';
%     [ISource, ITarget, ISourceNotInTarget] = ISourceTo(C1, C2)
%     C1(ISource) %duplicate free data in target order
%     C2(ITarget) %available data from target set    
%     C1(ISourceNotInTarget) %items without representative in target array

%     clc
%     C1 = {DB.Signatures.NFkappaB.probesMetaInfo.UniGene}';
%     C2 = {DB.probesMetaInfo.UniGene}';
%     [ISource, ITarget, ISourceNotInTarget] = ISourceTo(C1, C2)
%     C1(ISource) %duplicate free data in target order
%     C2(ITarget) %available data from target set    
%     C1(ISourceNotInTarget) %items without representative in target array
