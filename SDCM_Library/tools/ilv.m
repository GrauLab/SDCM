%ABSTRACT
%  Inline replacement for variables; for assignments in an inline Matlab expression.
%SYNTAX 
%  [result1, result2, ...] = ilv(variable1, variable2, ..., fcns, pipeFcns1, pipeFcns2, ...)
%INPUT
%  variables: any Matlab expression.
%  fcns: a single function handle or a cell array of function handles with the signature
%        @(variable1, variable2, ...)functionBody(...)
%  pipeFcns1: if provided, after calculating [results{:}]=fcns(variables{:}), another
%        ilv command is piped as: [pipeResult1, pipeResult2, ...] = ilv(result1, result2, ...
%        pipeFcns1, pipeFcns2, ...), i.e. this effectively allows writing a Matlab program in a single line.
%        Note that in this case the nargin(pipeFcns1) determines the requested nargout(fcns)
%        and the nargout of this ilv call will be the number of outputs requested from the last pipeFcnsJ.
%OUTPUT
% [result1, result2, ...]: the results requested from fcns(variables{:});
%EXAMPLES
%  In Matlab you cannot write expressions like "@(a)temp=longExpression(a);temp^2-log(temp)".
%  With inlineVar you can write instead: @(a)ilv(longExpression(a),@(x)x^2-log(x))
%  You can even write small programs in a single line with the pipeFcns options:
%     expectation = ilv(...
%       logspace(-10,0,10000) ... %set X axis.
%      ,@(X)ilsub({normpdf(X,0.5,0.1).*X, (diff(X([1,1:end]))+diff(X([1:end,end])))/2},':',{}) ... %calculate integrands and DX
%      ,@(integrands,DX)sum(integrands.*DX) ... %numeric integration
%     )
%  Real world examples:
%   a)  @(beadResults) ilv(...
%         inlineGetNthOutput([2,3,4], ttest2(vertcat(beadResults{IPheno1}), vertcat(beadResults{IPheno2}), 0.01, 'both', 'unequal')) ...
%        ,@(results)struct('p',results{1},'ci',results{2},'sd',results{3}.sd,'t',results{3}.tstat) ...
%       )...
%   b) commonQuantileValues_bgCorrected = parcellfun(...
%        slimFcn(@(x,xi)ilv(...
%          normpdf(x,S+mu,sigma) .* PDF_Signal_startingAtZero / PDF_Gemessen(xi) ... %=PDF_SignalWennGemessen_forx
%         ,@(PDF_SignalWennGemessen_forx) PDF_SignalWennGemessen_forx / sum(PDF_SignalWennGemessen_forx .* DS) ... %Normierung
%         ,@(PDF_SignalWennGemessen_forx_normed) sum(PDF_SignalWennGemessen_forx_normed.*S .* DS) ... %numerische Integration 
%        ),'mu','sigma','X','S','DS','PDF_Signal_startingAtZero','PDF_Gemessen') ... %vars from this workspace required on worker side.
%       ,num2cell(commonQuantileValues) ...
%       ,num2cell(round(interp1(X,1:length(X),commonQuantileValues,'nearest','extrap'))) ... %xi into PDF_Gemessen (table lookup)
%       ,true ... %combine results.
%      );
%NOTE:
%- Use slimFcn and specify the workspace variables needed by the function handle explicitly whenever
%  you plan to execute your command on workers. Else the process is inefficient or you might even run out of memory.
%  To disambigue the calling syntax, all input arguments starting from the last that are function handles 
%  (or cell arrays thereof) will be associated with ..., fcns, pipeFcns1, ... I.e. in order to pass a function handle
%  (or a cell array of function handles) as a variable_i, it must be followed by a non-function_handle variable,
%  for example a dummy variable [] as in: ...=ilv(myVar1,@(a)a^2,[],fcns,...).
%  Note that ilv checks that the count of input variables and nargin(fcns) match and errors, if it does not.
%- If you want to pass a function handle as variable you must wrap it as {{fcn}} and unwarp it from within the function that should be executed on it.
%AUTHOR
%  (C) Michael Grau, 2012
%  Library function for SDCM standalone deployment. Part of the Matlab inline 
%  language toolbox. (Not to be redistributed separately.)

function varargout = ilv(varargin)
  if(nargin<1) error('Not enough input arguments; see help.'); end
  %initialize:
    %parse inputs: get the ..., fcns, pipeFcns1, pipeFcns2, ... input arguments:
      BIsFcnsInput = cellfun(@(c)~isempty(c)&&(isa(c,'function_handle')||iscell(c)&&all(ilsub(cellfun(@(c)isa(c,'function_handle'),c),':'))), varargin);
      IIsFcnsInput = ilsub(find(BIsFcnsInput),@(self)self>=max([0,find(~BIsFcnsInput,1,'last')])+1);
      assert(length(IIsFcnsInput)>=1, 'missing input fcns and optional pipeFcns1, ...; see help.');
      fcns = varargin{IIsFcnsInput(1)};
      variables = varargin(1:IIsFcnsInput(1)-1); 
    %signature check:
      if(iscell(fcns))
        nExpectedInputs = unique(cellfun(@nargin, fcns));
        assert(isscalar(nExpectedInputs), 'If fcns is a cell array of function handles, all function handles must have identical nargin.');
      elseif(isa(fcns,'function_handle'))
        nExpectedInputs = nargin(fcns);
      else
        error('fcns must be a function handle or a cell array of function handles.');
      end
      if(nExpectedInputs>=0 && nExpectedInputs~=length(variables)) %NOTE: nExpectedInputs<0 is returned for varargin functions.
        error(...
           'Function "%s" expected %d variables, but %d were input to ilv. (If you need function handles as variables, make sure you add a dummy variable before the function input arguments to ilv; see help.)'...
          ,iif(isa(fcns,'function_handle'), @()mapSubPath(functions(fcns),'function'), sprintf('(cell array of %d function handles)',numel(fcns))) ...
          ,nExpectedInputs ...
          ,length(variables) ...
        );
      end
  %Apply fcns to variables and get all results:
    if(length(IIsFcnsInput)==1) %the results from fcns shall be output by this ilv call
      varargout = cell(1,nargout);
    else %piping case: the next function handle in the pipe determines the number of requested outputs from fcns:
      if(isa(varargin{IIsFcnsInput(2)},'function_handle'))
        nargin4NextFcn = nargin(varargin{IIsFcnsInput(2)});
        if(nargin4NextFcn<0)
          warning(...
             'The number of input arguments for the next function in the pipe %s could not be determined and is possibly variable.\nThus the number of output arguments to request from the predecessor function %s in the pipe is unclear.\nUse explicit @(x,y)fcnWithVariableInputs(x,y) syntax to disambigue this case.\nNow defaulting to one requested output.'...
            ,strtrim(strrep(evalc('varargin{IIsFcnsInput(2)}'),sprintf('\n'),''))...
            ,strtrim(strrep(evalc('fcns'),sprintf('\n'),''))...
          );
          nargin4NextFcn = 1;
        end
        varargout = cell(1,nargin4NextFcn); 
      else %in this case it must be a cell array of functions:
        nargin4NextFcn = nargin(varargin{IIsFcnsInput(2)}{1}); 
        if(nargin4NextFcn<0)
          warning(...
             'The number of input arguments for the next function in the pipe %s could not be determined and is possibly variable.\nThus the number of output arguments to request from the predecessor function %s in the pipe is unclear.\nUse explicit @(x,y)fcnWithVariableInputs(x,y) syntax to disambigue this case.\nNow defaulting to one.' ...
            ,strtrim(strrep(evalc('varargin{IIsFcnsInput(2)}{1}'),sprintf('\n'),''))...
            ,strtrim(strrep(evalc('fcns'),sprintf('\n'),''))...
          );
          nargin4NextFcn = 1;
        end
        varargout = cell(1,nargin4NextFcn); 
      end
    end
    if(iscell(fcns))
      if(~isempty(varargout))
        [varargout{:}] = cellfun(@(fcn)fcn(variables{:}), fcns, 'UniformOutput',false);
      else
        cellfun(@(fcn)fcn(variables{:}), fcns, 'UniformOutput',false); %allow command syntax ilv(...,@command) without any output arguments.
      end
    elseif(isa(fcns,'function_handle'))
      try
        if(~isempty(varargout))
          [varargout{:}] = fcns(variables{:});
        else
          %%warning('no output arguments requested from ilv; will run function "%s" in command mode returning zero outputs.',strtrim(evalc('disp(fcns);')));
          %fcns(variables{:}); %allow command syntax ilv(...,@command) without any output arguments.
          result4ans = fcns(variables{:}); 
          varargout = {result4ans};  
        end
      catch ex
        warning('the error reported below occurred in user function\n    ''%s''\nfor input %s', ilsub(functions(fcns),'.function'), evalc('variables'));
        rethrow(ex);
      end
    end
  %pipe results to additional pipeFcns recursively, if provided:
    if(length(IIsFcnsInput)>1)
      midPipeResults = varargout;
      if(nargout>0)
        varargout = cell(1,nargout);
      else
        varargout = cell(1,1); 
      end
      [varargout{:}] = ilv(midPipeResults{:}, varargin{IIsFcnsInput(2:end)});
    end
end
