%ABSTRACT:
%  Inline switch.
%SYNTAX: 
%  ils(actualValue, compareValue1, result1, compareValue2, result2, ...)
%  ils(actualValue, compareValue1, result1, compareValue2, result2, ..., customCompareFcn)
%NOTES: 
% - the last compareValue can be 'otherwise' to specify a fallback value
% - all (compareValue,result) pairs after the compare value 'otherwise' are ignored.
% - if the input agrs have even count and the last input arg is a customCompareFcn,
%   for example @(actualValue,candidateValue)~isempty(strfind(actualValue,candidateValue))
%   then this compare function is used instead of the standard.
% - every result_i value can be a @() function handle to only evaluate the result_i that will be returned from ils.
% - it can also be a @(actualValue) function handle, which is especially useful for the otherwise clause.
% - if actualValue is a singleton cell array, ils is mapped over the first non-singleton dimension 
%   (to apply ils to a cell array without mapping wrap it with additional {}).
%   The compareValues and resultValues are passed on to the mapped ils functions unchanged (there is no mapping of the compareValues any longer due to an ambiguity problem.)
%EXAMPLES
%  ils('zwei','eins',1,'zwei',2,'drei',3,'otherwise',nan)
%  ils({'zwei','eins','neun','drei'},'eins',1,'zwei',2,'drei',3,'otherwise',@()ilx(@warning,'otherwise clause invoked'))
%  ils('my string containing a substring','taining',true,'not',false,@(s,p)~isempty(strfind(s,p)))
%AUTHOR
%  (C) Michael Grau, 2012
%  Library function for SDCM standalone deployment. Part of the Matlab inline 
%  language toolbox. (Not to be redistributed separately.)

function varargout = ils(actualValue, varargin)
  bIsMappingRecursion = ischar(varargin{end}) && strcmp(varargin{end},'isRecursion');
    if(bIsMappingRecursion)
      varargin = varargin(1:end-1);
    end
    if(~bIsMappingRecursion && iscell(actualValue))
      iFirstNonSingleTonDim = find(size(actualValue)>1,1);
        if(isempty(iFirstNonSingleTonDim)) iFirstNonSingleTonDim = 1; end
      if(sum(size(actualValue)>1)>1)
        actualValue = num2cell(actualValue,setdiff(1:length(size(actualValue)),iFirstNonSingleTonDim));
      end
      bSingleThreaded = numel(actualValue)<500;      
      if(~bSingleThreaded)
        %NOTE: only works for scalar resultValue: result = parcellfun(slimFcn(@(value)ils(value, varargin{:}, 'isRecursion'), 'varargin'), actualValue);
        %expand the compareValue strings and scalar result values to the input array size (in oder to support non-scalar resultValueArrays):
          for i=1:length(varargin)
              varargin{i} = repmat(varargin(i),size(actualValue));
          end
          result = parcellfun(slimFcn(@(value, varargin) ils(value, varargin{:}, 'isRecursion')), actualValue, varargin{:}, false, true, bSingleThreaded);
      else %faster loop variant:
        result = cell(size(actualValue));
        for i=1:numel(actualValue)
          result{i} = ils(actualValue{i}, varargin{:}, 'isRecursion');
        end
      end
      varargout = {result};
      return;
    end
  bHasCustomFcnComparer = mod(length(varargin),2)==1 && isa(varargin{end},'function_handle') && nargin(varargin{end})==2;
    if(bHasCustomFcnComparer)
      customCompareFcn = varargin{end};
      varargin = varargin(1:end-1);
    else
      if(isa(varargin{end},'function_handle') && nargin(varargin{end})==2 && ~any(cellfun(@(resultValue)isa(resultValue,'function_handle'),varargin(2:2:end-1))))
        warning('the last input appears to be a custom compare function, but is used as resultValue for the (end-1)-input case argument; ist this is intended to be a custom compare function, check that the total number of input arguments is even.');
      end
      customCompareFcn = [];
    end
  if(mod(numel(varargin),2)~=0) error('inline switch function syntax error; see help'); end; 
  iOtherwiseClause = (find(strcmp('otherwise',varargin(1:2:end)), 1, 'last')-1)*2+1;
  bHasOtherwiseClause = ~isempty(iOtherwiseClause);
    if(bHasOtherwiseClause)
      otherwiseResult = varargin{iOtherwiseClause+1};
    end
    
  for i=1:numel(varargin)/2
    compareValue = varargin{(i-1)*2+1};
    if((i-1)*2+1==iOtherwiseClause) continue; end; %otherwise will be handles below.
    bTreatErrorAsFalse = false; %older behavior was bHasOtherwiseClause, but even with otherwise clause, comparison errors should not be silent.
    if(isEqual(actualValue, compareValue, bTreatErrorAsFalse, customCompareFcn))
      result = varargin{(i-1)*2+2};
      if(isa(result,'function_handle'))
        if(nargin(result)==0)
          %result = result();
          if(nargout==0)
            warning('ils was called and the return value has to be determined by a function call, but zero outputs are requested; the function handle for the return vale will be called now, but its output (if any) will not be requested nor returned. Consider adding a left-hand side variable in the call like myResult=ils(...).');
            result();
          else
            varargout = cell(1,nargout);
            [varargout{:}] = result();
          end
        elseif(nargin(result)==1)
          %result = result(actualValue);
          if(nargout==0)
            warning('ils was called and the return value has to be determined by a function call, but zero outputs are requested; the function handle for the return vale will be called now, but its output (if any) will not be requested nor returned. Consider adding a left-hand side variable in the call like myResult=ils(...).');
            result(actualValue);
          else
            varargout = cell(1,nargout);
            [varargout{:}] = result(actualValue);
          end
        else
          error('resultValue "%s" was a function handle of unsupported signature; if the function handle itself should be the result, prefix a @().', evalc('result'));
        end
      else
        varargout = {result};
      end
      return;
    end
  end
  if(bHasOtherwiseClause)
    result = otherwiseResult;
    if(isa(result,'function_handle'))
      if(nargin(result)==0)
        %result = result();
        varargout = cell(1,nargout);
        [varargout{:}] = result();
      elseif(nargin(result)==1)
        %result = result(actualValue);
        varargout = cell(1,nargout);
        try
          [varargout{:}] = result(actualValue);
        catch ex
          switch(ex.identifier)
            case 'MATLAB:assigningResultsIntoInitializedEmptyLHS'
              error('ils was called with zero output arguments, but user results function "%s" returned more', ilsub(functions(result),'.function'));
            otherwise
              rethrow(ex);
          end
        end
      else
        error('resultValue "%s" was a function handle of unsupported signature; if the function handle itself should be the result, prefix a @().', evalc('result'));
      end
    else
      varargout = {result};
    end
  else %error, if there is no otheriwse branch, because in this scenario, the result is undefined:
    error('{%s} was not equal any of the provided compareValues={%s} and there was no ''otherwise'' value specified.'...
      ,iif(ischar(actualValue),actualValue,evalctrim('actualValue')) ...
      ,evalctrim('varargin(1:2:end)') ...
    ); 
  end
end

function b = isEqual(actualValue, compareValue, bTreatErrorAsFalse, customCompareFcn)
  try
    if(~isempty(customCompareFcn))
      b = customCompareFcn(actualValue, compareValue);
    elseif(isa(compareValue,'function_handle')&&nargin(compareValue)==1) %allow using predicate functions as compare value, e.g. ils(IC50, @(x)x<0.01,@(x)num2str(x,'%0.2e'), @isinf,'(insensitive)', 'otherwise',@(x)num2str(x,'%0.2f'))
      try
        b = compareValue(actualValue);
      catch ex
        ex.addCause(sprintf('Custom compare predicate "%s" caused an error for actual value = %s.', ilsub(functions(compareValue),'.function'), strtrim(evalc('actualValue'))));
        rethrow(ex);
      end
    elseif(isnumeric(actualValue)||islogical(actualValue))
      if(~(isnumeric(compareValue)||islogical(compareValue))) error([evalc('actualValue'), ' is numeric, but ',evalc('compareValue'),' is not => cannot compare.']); end
      if(~all(size(actualValue)==size(compareValue))) error(['size of ', evalc('actualValue'), ' does not match size of ',evalc('compareValue'),' => cannot compare.']); end
      b = all( (actualValue(:)==compareValue(:)) | (isnan(actualValue(:))&isnan(compareValue(:))) );
    elseif(ischar(actualValue))
      if(~ischar(compareValue)) error([evalc('actualValue'), ' is a string, but ',evalc('compareValue'),' is not => cannot compare.']); end
      b = strcmp(actualValue, compareValue);
    else
      error(['comparison of ', evalc('actualValue'), ' with ',evalc('compareValue'),' is not implemented.']);
    end
  catch ex
    if(bTreatErrorAsFalse)
      b = false; return;
    else
      ex.rethrow();
    end      
  end
end

function sResult = evalctrim(sEval)
  sResult = evalin('caller', ['evalc(''',sEval,''')']);
  sResult = strrep(sResult,'ans =','');
  sResult(sResult==10 | sResult==13) = [];
  l = length(sResult); bShortened = true;
  while(bShortened)
    sResult = strrep(sResult,'  ', ' ');
    bShortened = length(sResult)<l;
    l = length(sResult);
  end
end
