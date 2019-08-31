%ABSTRACT
% Inline if for anonymous functions.
%SYNTAX
% result = iif(condition, trueResult, falseResult)
%INPUT
% condition: a single boolean scalar. Can be a boolean array; in this case
%            scalar results get repmatted to the same size; if the results are
%            scalar neither, they need to be of the same size as condition.
%            Can also be a function handle of signature @(trueResult,falseResult)condition 
%            or @(trueResult)condition, i.e. the condition gets calculated based on 
%            trueResult and falseResult. (cf. NOTE2)
% trueResult: any Matlab value or a function handle with no arguments
%             returning any Matlab value (single output; see NOTE1)
% falseResult: any Matlab value or a function handle with no arguments
%             returning any Matlab value (single output; see NOTE1)
%OUTPUT
% result: trueResult, if condition, elese falseResult.
%NOTE1
% Sometimes you do not want that e.g. the falseResult expression gets calculated,
% if the condition is true. Usually Matlab calculates everything that is on the
% line and then calls the iif function. To prevent this, you simply need
% to prefix the expression by a "@()", i.e. you wrap it as an anonymous function
% that takes no arguments and returns a single output.
% EXAMPLES
%   iif(true, 1, 'hello') yields 1.
%   iif(false, 1, 'hello') yields 'hello'.
%   iif([true,false], [1,2], [3,4]) yields [1,4].
%   iif([true;false], {'a';'b'}, {'c';'d'}) yields {'a';'d'}.
%   iif(true, 1, ilerr('test')) causes an error to be thrown, although it is in the false branch.
%   iif(true, @()1, @()ilerr('test')) only calculates the true branch and returns 1.
%   iif(iscell(c)&&length(c)==1,@()c{1},c): if c is e.g. {1}, return 1, if it is {} return as is without error.
% Note that only function handles with nargin(fncHandle)==0 (i.e. "@()..." type of handles) will be evaluated,
% i.e. iif(true, @(a)b, @(c)d) will return @(a)b and not try to call this handle without any input arguments.
% This allows you to use iif also with function handles as results. If one iif branch needs to return
% a function handle with zero input arguments, you must precede it by an additional "@()":
%   iif(true, @()@()temp, false): returns the @()temp function handle.
%NOTE2
% Sometimes you need to stance out certain elements, i.e. have a locial array as condition
% that is based on trueResult or falseResult. In this case you can provide a function handle (a predicate in logical terms)
% as the condition input parameter. This syntax is especially convenient, if the term for trueResult and/or falseResult are long-in-code.
% EXAMPLE
%   iif(@(tR)~isnan(tR), tR, 0) %returns an array of size tR with all NaNs replaced by zeros.
%   iif(@(tR,fR)tR>0 & fR<0, tR, fR) %returns tR elements, if positive and the matiching fR element it negative; else returns fR elements instead.
% Note that when tR or fR are given in function handle format (see NOTE1) your
% condition function handle will receive their evaluated form. As providing trueResult or falseResult
% as function handle has the sole purpose of preventing them from being evaluated before the condition has been evaluated,
% they could as well be provided directly in this scenario (as they are required to calculate the conditions).
%AUTHOR
% (C) Michael Grau, 2012
% Library function for SDCM standalone deployment. Part of the Matlab inline 
% language toolbox. (Not to be redistributed separately.)

function result = iif(condition, trueResult, falseResult)
  if(nargin<3) %useful for inline version of "if(condition)fncBody();end"
    if(isnumeric(trueResult) && isscalar(trueResult))
      falseResult = NaN; 
    elseif(ischar(trueResult))
      falseResult = ''; 
    else
      error('no default falseResult exists for this trueResult');
    end
  end 
  if(islogical(condition))
    if(length(condition)==1) %scalar condition:
      if(condition)
        result = trueResult;
      else
        result = falseResult;
      end
      if(isa(result,'function_handle') && nargin(result)==0)
        result = result();
      end
    else %vektoriell:
      %calculate needed result branche(s) provided as function handles:
        bTrueBranchNeeded = any(condition(:));
        bFalseBranchNeeded = any(~condition(:));
        if(isa(trueResult,'function_handle') && nargin(trueResult)==0 && bTrueBranchNeeded)
          trueResult = trueResult();
        end
        if(isa(falseResult,'function_handle') && nargin(falseResult)==0 && bFalseBranchNeeded)
          falseResult = falseResult();
        end
      %if the results are scalar, repmat them to condition size:
        if(length(trueResult)==1) trueResult=repmat(trueResult,size(condition)); end;
        if(length(falseResult)==1) falseResult=repmat(falseResult,size(condition)); end;
      %assert compatible sizes:
        assert(all(size(condition)==size(trueResult)), 'non-scalar condition and non-scalar trueResult with different array sizes are invalid syntax for iif');
        assert(all(size(condition)==size(falseResult)), 'non-scalar condition and non-scalar falseResult with different array sizes are invalid syntax for iif');
      %return result using condition as boolean pattern/stamp matrix:
        if(isnumeric(trueResult) && isnumeric(falseResult))
          result = nan(size(condition));
            result(condition) = trueResult(condition);
            result(~condition) = falseResult(~condition);
        elseif(iscell(trueResult) && iscell(falseResult))
          result = cell(size(condition));
            result(condition) = trueResult(condition);
            result(~condition) = falseResult(~condition);
        else
          error('non-scalar condition requires trueResult and falseResult to both be either numeric arrays or cell arrays');
        end
    end
  elseif(isa(condition,'function_handle'))
    if(nargin(condition)==1)
      if(isa(trueResult,'function_handle') && nargin(trueResult)==0)
        trueResult = trueResult();
      end
      condition = condition(trueResult);
    elseif(nargin(condition)==2)
      if(isa(trueResult,'function_handle') && nargin(trueResult)==0)
        trueResult = trueResult();
      end
      if(isa(falseResult,'function_handle') && nargin(falseResult)==0)
        falseResult = falseResult();
      end
      condition = condition(trueResult, falseResult);
    end
    result = iif(condition, trueResult, falseResult);
  else
    error('unsupported type for input parameter "condition"');
  end
end
