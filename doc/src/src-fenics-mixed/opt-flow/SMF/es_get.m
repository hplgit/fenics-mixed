function data = es_get(es, strname)
% function data = es_get(es [, strname])
%   es_get provides data from es, e.g. the result of an optimization 
%
% Input arguments:
%   es : evolution strategy struct, see function es_new
%
%   [strname] : (optional, default is 'result') Can be one of the 
%               following.
%
%           [] : If strname is empty, the default 'result' is used
%                (see below). 
% 'alloptions' : data is an array of structs with fields iter
%                and options. The field iter holds the iteration
%                number, when (new) options were set. The field
%                options is an options struct (see es_options_new),
%                which holds the control parameters for the evolution
%                strategy.
%      'input' : data is a struct with fields N, xlow and xup, which hold
%                the input arguments from calling es_new, and field
%                options. This is all the information needed to
%                replicate previous results. Example:
%                   es1 = es_new(N, xlow, xup);
%                   stinput = es_get(es1, 'input');
%                The following three calls are mutually identical
%                   es2 = es_new(es1, 'clone');
%                   es2 = es_new(N, xlow, xup, es_get(es1, 'options'));
%                   es2 = es_new(stinput);
%                The following two calls are mutually identical
%                   es3 = es_new(es1);
%                   es3 = es_new(N, xlow, xup, es_options_new('es', es1));
%                es1 and es2 are completely identical, while es3
%                runs with different random number realizations. 
%    'options' : data is an options struct (see es_options_new),
%                which holds the control parameters for
%                the evolution strategy. 
%     'result' : (default) data is a structure with fields
%                  funevals : Number of executed function evaluations
%                  funvalue : Function value of best search point
%                             of current iteration. 
%                  x : Best evaluated search point of current iteration.
%                  xreco : Weighted recombination of best search
%                          points of current iteration. Unlike x,
%                          xreco may lie outside of the given
%                          bounds. Even thought not evaluated, xreco
%                          is expected to have a better function value
%                          than x.
%                  best : struct with fields
%                    funevals : Number of executed function
%                               evaluations when the best evaluated
%                               search point best.x occurred.
%                    funvalue : Function value of search point best.x
%                    x : Best evaluated search point found so far. 
%                  tolx : coordinate wise mean change (vector 
%                         es.sigma*sqrt(diag(es.C))). 
%                  axes : axis lengths of the mutation ellipsoid
%                         disregarding step size sigma (diag(es.D)).
%                  sigma : step size es.sigma
%                  axisratio : ratio between longest and shortest
%                              axis of the mutation ellipsoid. 
%  'stopflags' : data is a sorted array of flag numbers consisting of
%                all active stop criteria that cause a break off (for
%                stop criteria see function es_options_new). If no
%                stop criterion is currently fulfilled, data is an
%                empty array, that is isempty(data) is true.
%          'x' : data is the current search point <x>. 
% 
% For examples call es_help (or see file es_READ.ME). 
%
% See also es_new, es_set 

%--------%--------%--------%--------%--------%--------%--------%--------%

  error(nargchk(1,2,nargin));
  % Set default
  if nargin == 1
    strname = 'result';
  elseif isempty(strname)
    strname = 'result';
  end
  
  % Check strname
  if strcmp('input', strname)
    data = es.input;
    data.options = es.aroptions(1).options;
  elseif strcmp('options', strname)
    data = es.options;
    if length(es.aroptions) > 1
      warning(['Options were changed during optimization, use' ...
	    ' es_get(es, ''alloptions'') to get complete information']);
    end
  elseif strcmp('alloptions', strname)
    data = es.aroptions;
  elseif strcmp('firstoptions', strname)
    data = es.aroptions(end).options;
  elseif strcmp('lastoptions', strname)
    data = es.options;
  elseif strcmp('result', strname)
    %@@@
    data.funevals = es.counteval;
    data.funvalue = es.rec.funval;
    data.x = es.rec.x;
    data.xreco = es.xreco;
    data.best.funevals = es.rec.best.eval;
    data.best.funvalue = es.rec.best.funval;
    data.best.x = es.rec.best.x;
    data.tolx = es.sigma * sqrt(diag(es.C));
    data.axes = diag(es.D);
    data.sigma = es.sigma;
    data.axisratio = max(data.axes)/min(data.axes);
    data.README = 'type "help es_get" for field details (view item ''result'')';
  elseif strcmp('stopflags', strname)
    data = sort(es.arflgstop);
  elseif strcmp('x', strname)
    data = es.xreco;
  else
    error([strname ' is not a valid input']);
  end




