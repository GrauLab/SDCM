%ABSTRACT
%  Allows you to call a function that returns multiple outputs inline and retrieves
%  only the output(s) as indicated by N.
%SYNTAX:
%  varargout = ilnth(N,fcn,varargin), with N={1,3,5} for multiple outputs or =[1,3,5] for a single 1*3 cell array wrapping output.
%INPUT
%  N: Specifies the output(s) of the fcn to return. The following possibilities are implemented:
%     - If N is an index cell vector {i1,i2,..}, inlineGetNthOutput returns the ith outputs like 
%       a function, i.e. the caller gets only the amout of outputs he requested.
%     - If N is a numeric index i: returns the ith output. (Identical to N={i} calling syntax.)
%       Note: To wrap a single output as a cell array of length 1, use {inlineGetNthOutput(i,fcn,..)}.
%     - If N is a numeric index vector [i1,i2,..], returns the ith outputs of fcn, but wraps
%       them in a cell array, i.e. only one output is returned from inlineGetNthOutput
%       with the value {i1th output, i2th output, ...}.
%  fcn: Function handle to the function to call.
%  varargin: Any parameters to pass on to fcn.
%AUTHOR
%  (C) Michael Grau, 2012
%  Library function for SDCM standalone deployment. Part of the Matlab inline 
%  language toolbox. (Not to be redistributed separately.)

function varargout = ilnth(N,fcn,varargin)
  %Default params:
    if(isnumeric(N)&&length(N)==1) N={N}; end; %return a single output unwrapped per default.
    assert(isa(fcn,'function_handle'),'fcn must be a function handle');
  %Initialize outputs array for fcn:
    if(iscell(N))
      maxN = max([N{:}]);
    else
      maxN = max(N);
    end
    outputs = cell(1,max(maxN));
  %call fcn:
    [outputs{:}] = fcn(varargin{:});
  %return outputs wrapped or as flat (multiple outputs):
    if(iscell(N))
      varargout = outputs([N{:}]); %cell-index into outputs, i.e. return multiple outputs of fcn as multiple outputs
    else
      varargout = {outputs(N)}; %array.index into outputs, i.e. return one output that is a cell array of selected outputs of fcn.
    end  
end

