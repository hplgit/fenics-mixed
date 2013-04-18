function es = es_set(es, strname, value)
% function es = es_set(es, strname, value)
%   es_set sets the parameter called strname to value
%
% Input Arguments: 
%   es : evolution strategy struct, see es_new. 
%
%   strname : Parameter name (string). Valid names are: 
%     'options' (type of value: options struct, see
%        es_options_new) : Sets new options in es. 
%     'sigmafac' (type of value: scalar) : Modifies overall
%        step-size by multiplication with value. 
%     'x' (type of value: N by one matrix): Sets current search
%        point <x> to value. The search continues from the new search
%        point.
%
% Example: 
%      es = es_new(10, -2, 1);
%      es = es_run(es, 'norm', 20);
%      es = es_set(es, 'x', ones(10,1)); % set new x-vector
%      es = es_run(es, 'norm', 20);
%      % get and modify the options
%      opts = es_get(es, 'options');
%      opts = es_options_set(opts, 'lbound', -1, 'ubound', 0.5);
%      % set the modified options
%      es = es_set(es, 'options', opts);                         
%
% See also es_new, es_get, es_options_new

%
% Remark: es_set(es, 'options') is called in es_new and therefore at
% least once for every es struct
%
%--------%--------%--------% BEGIN es_set %--------%--------%--------%--------%
  error(nargchk(1,3,nargin));
  if nargin == 2
    error(['value for "' strname '" not given']);
  elseif nargin == 3
    if strcmp('options', strname)
      es = es_setoptions(es, value);
    elseif strcmp('sigmafac', strname) 
      if value <= 0
	error('value for "sigmafac" must ge greater than zero');
      end
      es.sigma = es.sigma * value;
    elseif strcmp('x', strname)
      if ~isequal(size(value), [es.N 1])
	error('Wrong sized x-value');
      end
      es.xreco = value;  
    else
      error([strname ' is not a valid argument']);
    end
  end
  
  % Check consistency of es.opts and es.sp
  if ~isempty(find(es.opts.lbound >= es.opts.ubound)) ...
	| es.sp.cc > 1 | es.sp.cs >= 1 | es.sp.ccov >= 1 ...
	| es.sp.cc <= 0 | es.sp.cs <= 0 | es.sp.ccov < 0 ...
	| es.sp.damp <= 0 ...
	| es.sp.lambda < 2 ...
	| es.sp.lambda <= es.sp.mu ...
	| es.sp.lambda < length(es.sp.recoweights) ...
	| ~isequal(size(es.opts.xstart), [es.input.N 1]) ...
	| ~isequal(size(es.opts.scaling), [es.input.N 1]) ...
	| es.opts.updatemod < 0
    es.sp
    error(['Consistency check failed. Possibly ' ...
	  'mu>=lambda or ' ...
	  'lambda<length(recoweights) or ' ...
	  'other strategy parameter setting invalid or ' ...
	  'variable bounds incorrect or ' ...
	  '...???']);
  end
%--------%--------%--------% END es_set %--------%--------%--------%--------%

function es = es_setoptions(es, options)
%
% Interprets the options struct for the es struct. This is the
% interface between options and es. 
%
  
%--------%--------%--------% BEGIN es_setoptions %--------%--------%--------%
  N = es.input.N; % for readability only
  
  % Set stop criterions
  es.opts.maxfunevals = evalpara(options.maxfunevals, es);
  es.opts.tolx = evalpara(options.tolx, es);
  es.opts.tolfun = evalpara(options.tolfun, es);
  es.opts.minfunvalue = evalpara(options.minfunvalue, es);
  es.opts.flgnumprec = evalpara(options.numprecstop, es);
  es.opts.flgmaxcond = evalpara(options.maxcondstop, es);

  % Bounds and Minimal Step-Size
  if size(options.lbound, 2) ~= 1 | ...
	(size(options.lbound, 1) ~= 1 & ...
	size(options.lbound, 1) ~= N)
    error('Option lbound has unvalid size')
  end
  if size(options.ubound, 2) ~= 1 | ...
	(size(options.ubound, 1) ~= 1 & ...
	size(options.ubound, 1) ~= N)
    error('Option ubound has unvalid size')
  end
  es.opts.lbound = ones(N,1).*options.lbound;
  es.opts.ubound = ones(N,1).*options.ubound;
  es.opts.minxchange = ones(N,1).*options.minxchange; 

  % Strategy Parameters: Selection (lambda, mu, recoweights)
  lambda = evalpara(options.lambda, es);
  if lambda ~= es.sp.lambda
    es.arfunval = [];
    es.sp.lambda = lambda;
  end
  es.sp.mu = evalpara(options.mu, es);
  es.sp.recoweights = evalpara(options.recoweights, es);
  if size(es.sp.recoweights, 1) == 1 % recoweights is column-vector
    es.sp.recoweights = transpose(es.sp.recoweights);
  end
  if ~isempty(find(es.sp.recoweights <= 0))
    error('all recoweights must be greater than zero');
  elseif ~isequal(sort(es.sp.recoweights), ...
	es.sp.recoweights(end:-1:1)) 
    error('recoweights must be descending');
  elseif ~isequal(size(es.sp.recoweights), [es.sp.mu 1]) 
    error(['mu==' num2str(es.sp.mu) ' and length of recoweights (==' ...
	  num2str(length(es.sp.recoweights)) ') disagree']);
  end
  es.const.cw = sum(es.sp.recoweights)/norm(es.sp.recoweights);
  es.const.recow = es.sp.recoweights/sum(es.sp.recoweights);

  % Strategy parameters: Adaptation rates
  es.sp.cc = evalpara(options.cc, es);
  es.const.ccu = sqrt(es.sp.cc*(2-es.sp.cc));
  es.sp.cs = evalpara(options.cs, es);
  es.const.csu = sqrt(es.sp.cs*(2-es.sp.cs));
  es.sp.damp = evalpara(options.damp, es);
  es.sp.ccov = evalpara(options.ccov, es);

  % Miscellaneous 
  es.opts.updatemod = evalpara(options.updatemod, es); 
  
  % Handle initialization of es only in generation zero
  if es.countgen == 0
    xlow = ones(N, 1).*es.input.xlow; 
    xup = ones(N, 1).*es.input.xup;
  
    % Initial x-value
    if isempty(options.xstart) % default
      es.opts.xstart = ones(N,1).*((xup+xlow)/2);
    else
      if size(options.xstart) == [N 1]
	es.opts.xstart = options.xstart; 
      elseif size(options.xstart) == [1 N]
	es.opts.xstart = options.xstart'; 
      else
	error('xstart has wrong size or dimension');
      end
    end
    if sum(es.opts.xstart < xlow) | sum(es.opts.xstart > xup)
      warning('xstart not within the original initial search cubus');
    end
    es.xreco = es.opts.xstart;
    
    % Initial sigma and scaling
    es.opts.sigmafacstart = evalpara(options.sigmafacstart, es);
    if sum(size(es.opts.sigmafacstart) ~= [1 1])
      error('options.sigmafacstart must be a scalar');
    end
    es.sigma = min(0.5, 2*es.const.cw/N^0.5)*es.opts.sigmafacstart;
    
    if isempty(options.scaling) % default
      es.opts.scaling = (xup-xlow)/2; % Nx1-vector
    else
      if size(options.scaling) == [N 1]
	es.opts.scaling = options.scaling;
      elseif size(options.scaling) == [1 N]
	es.opts.scaling = options.scaling';
      else
	error('scaling has wrong size or dimension');
      end
    end
    es.D = diag(es.opts.scaling);
    es.BD = es.B*es.D; 
    es.C = es.BD*transpose(es.BD);
    if cond(es.D) > 1e5
      warning('Initial search cubus is badly scaled');
      error(['To force the scaling suppress this error by editing function' ...
	    ' es_set']);
    end

    % Set random seed
    es.opts.randomseed = evalpara(options.randomseed, es);
    randn('state', es.opts.randomseed);
    
    % Miscellaneous 
    es.arflgstop = [];
    es.opts.initphase = options.initphase;
    if ~ismember(es.opts.initphase, [0 1]) 
      error('Option ''InitPhase'' must be 0 or 1');
    end
    es.flginitphase = es.opts.initphase; 
  end    

  % Save new options
  es.options = options;
  if ~isfield(es, 'aroptions') % i.e. es.countgen == 0
    es.aroptions(1).iter = es.countgen;
    es.aroptions(1).options = options; 
  else
    if es.aroptions(end).iter < es.countgen
      es.aroptions(end+1).iter = es.countgen;
    end
    es.aroptions(end).options = options;
  end

%--------%--------%--------% END es_setoptions %--------%--------%--------%

function res = evalpara(value, es)
%--------%--------%--------%--------%--------%--------%--------%--------%
  N=es.input.N; 
  if isfield(es, 'sp')
    if isfield(es.sp, 'lambda')
      lambda=es.sp.lambda; 
      if lambda > 1 & isfield(es, 'opts')
	if isfield(es.opts, 'maxfunevals')
	  maxiter = es.opts.maxfunevals / lambda;
	end
      end
    end
    if isfield(es.sp, 'mu')
      mu=es.sp.mu;
    end
    if isfield(es.sp, 'cs')
      cs = es.sp.cs;
    end
    if isfield(es.sp, 'ccov')
      ccov = es.sp.ccov;
    end
  end
  if ischar(value)
    res = eval(value);
  else
    res = value;
  end
  
  if isempty(res)
    error(['value ''' value ''' invalid']);
  end
%--------%--------%--------%--------%--------%--------%--------%--------%

