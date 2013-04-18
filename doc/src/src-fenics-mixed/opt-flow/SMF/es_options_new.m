function options = es_options_new(varargin)
% function options = es_options_new([varargin])
%   generates an options structure
% 
% The options struct is the user interface to the evolution strategy
% struct (see es_new).  es_options_new returns a struct with all
% (default) options for the evolution strategy.
%
% Input Arguments: 
%   varargin : (optional) Like ('param1', value1, 'param2', value2,...)
%              pairs of parameter names, which are fieldnames of the
%              options struct (additionally 'es' and 'options' are
%              valid, compare 6) below), and values which modify the
%              default value. The fieldname must be a string, the
%              value can be the corresponding value or can be a string
%              which will be evaluated with EVAL. Use upper case N to
%              indicate problem dimension. See examples below. The
%              following parameter names are valid (lower/upper case
%              does not matter):
%
%    1) Stop Criterions
%
%    'MaxFunEvals' (type of value: number, default==inf, where N
%      is the problem dimension): Stop criterion No.1. Maximal
%      number of function evaluations to be executed.
%
%    'TolX' (type of value: scalar or N by 1 matrix, default==1e-7): Stop
%      criterion No.2. If variation of search points becomes
%      considerably smaller than TolX in all coordinates, search is
%      stopped. TolX can only take effect if MinXChange is smaller
%      than TolX.
%
%    'TolFun' (type of value: scalar, default==0): Stop criterion
%      No.3. If variation of function value is smaller TolFun, search
%      is stopped.
%
%    'MinFunValue' (type of value: scalar, default==-inf): Stop
%      criterion No.4. If the function value of an evaluated search
%      point is below MinFunValue search is stopped. 
%
%    'NumPrecStop' (type of value: flag with value -3, -2, -1, 0, or 1;
%      default==0): Stop criterion No.5. If the algorithm detects a
%      numeric precision problem depending on the value of
%      'NumPrecStop' the following happens:
%      'NumPrecStop' == 1 : A warning is printed and the search is
%          stopped. 
%      'NumPrecStop' == 0 : (default) A warning is printed and the
%          step size is enlarged.
%      'NumPrecStop' == -1 : The step size is enlarged, no warning
%          is printed.
%      'NumPrecStop' == -2 : The problem is ignored, but a warning is 
%          printed.
%      'NumPrecStop' == -3 : The problem is ignored, no warning is 
%          printed.
% 
%    'MaxCondStop' (type of value: flag with value -3, -2, -1, 0, or 1;
%      default==0): Stop criterion No.6. If the maximal condition
%      number of the covariance matrix of 1e14 is exceeded depending
%      on the value of 'MaxCondStop' the following happens:
%       'MaxCondStop' == 1 : A warning is printed and the search is
%          stopped. 
%       'MaxCondStop' == 0 : (default) The condition number is 
%          is limited to 1e14 + 1 and a warning is printed. 
%       'MaxCondStop' == -1 : The condition number is 
%          is limited to 1e14 + 1, no warning is printed. 
%       'MaxCondStop' == -2 : Only a warning is printed. 
%       'MaxCondStop' == -3 : (not recommended) nothing is done.
%
%    2) Initial Values and Scaling
%
%    'XStart' (type of value: N by 1 matrix): Initial search
%      point. If no initial point is given, the search starts in the
%      middle of the given initial search hyper-cube (see
%      es_new). Remark that the initial search point is not evaluated
%      by the search strategy.
%
%    'SigmaFacStart' (type of value: positive scalar, default==1): The
%      initial step-size is multiplied by SigmaFacStart. That is, the
%      initial step lengths can be enlarged with values greater than
%      one and reduced with values smaller than one. For a given es
%      struct the complete initial distribution of the steps can be
%      determined by the options SigmaFacStart and Scaling.
%
%    'Scaling' (type of value: N by 1 matrix): Initial scaling of the
%      coordinate axes. Default is ones(N,1).*(xup-xlow)/2, where xup
%      and xlow are the mandatory parameters as given with es_new (see
%      there). This determines the initial step lengths of the
%      search. Changing 'Scaling' changes the (initial) step lengths.
%
%    'RandomSeed' (type of value: number or string which can be evaluated
%      with function eval): Initial value for pseudo random number
%      generator. If no value is given, the initial value is different
%      with each call of es_options_new.
%
%    3) Bounds and Minimal Step-Size
%
%    'LBound' (type of value: scalar or N by 1 matrix): Lower bound
%      for search. The strategy evaluates only points X with
%      X>=LBound. Default is -inf (no bound). 
%
%    'UBound' (type of value: scalar or N by 1 matrix): Upper bound for
%      search. Only points X with X<=UBound are evaluated.  Default is
%      inf (no bound).
%
%    'MinXChange' (type of value: scalar or N by 1 matrix,
%      default==0): Smallest sensible variation in each coordinate
%      axis. This minimal variation is ensured by adjustments of the
%      step size and the covariance matrix. Remark that MinXChange can
%      influence the impact of stop criteria TolX and TolFun.
%
%    4) Strategy Parameters
%
%    'lambda' : Number of offspring, default is '4+floor(3*log(N))':
%      N      = 2    3    4    5    6    7    8   10   12   20   50  100
%      lambda = 6    7    8    8    9    9   10   10   11   12   15   17
%      Setting lambda you can refer to N as in the default. Refering
%      e.g. to mu can lead to errors difficult to trace back. 
%      Remark that the defaults of mu and recoweights depend on
%      lambda. If only lambda was set, this alway leads to reasonable
%      values for mu and recoweights.
%
%    'mu': Number of parents, default is 'floor(lambda/2)'. Setting mu
%      you can refer to lambda as in the default or to N. If mu and no
%      recoweights are given, all recoweights are set to 1/mu which
%      implements the typical (mu,lambda)-ES with intermediate multi-
%      recombination (which is not default). It must hold mu<lambda.
%
%    'recoweights' (type of value: mu by 1 or 1 by mu matrix):
%      Weights for weighted recombination, default is
%      'log((lambda+1)/2)-log(1:mu)'. Setting recoweights you can
%      refer to mu and/or lambda as in the default and to N. It must
%      hold length(recoweights)<=lambda and mu==length(recoweights).
%
%    Furthermore the following identifyers can be used (defaults in
%    parentheses): 'cc' ('4/(N+4)'), 'ccov' ('2/(N+2^0.5)^2'), 'cs'
%    ('4/(N+4)'), 'damp' ('(1-min(0.7,N/maxiter))/cs + 1', where
%    maxiter equals MaxFunEvals/lambda)
%
%    5) Miscellaneous 
%
%    'InitPhase' (type of value: flag, default==1): By default an initial
%      phase is implemented where only global step-size adaptation
%      takes place. Adaptation of the covariance matrix begins after
%      a number of iterations if step-size appears to be sufficiently
%      large. This tries to circumvent problems connected with a
%      small initial step-size. If initphase is set to 0 adaptation
%      starts with the first iteration. 
%
%    'UpdateMod' (type of value: positive number, default==1):
%      Number of iterations done without update of the mutation
%      distribution. If strategy internal operations take a
%      considerable time compared to the time taken by objective
%      function evaluations it may help to set updatemod to a larger
%      value. In this case '1/ccov/N/5' (or more easily 'N/10' or
%      'sqrt(N)') is suggested. If updatemod==0 no update takes place
%      at all. (The covariance matrix is updated if mod(iteration,
%      updatemod) < 1).
%
%    6) An options struct or an es struct. This must be the first
%       argument in the list.
%
%       'options' (type of value: options struct): Options from the
%         input options struct are used with a new initial value for the
%         pseudo random number generator.
%
%       'es' (type of value: es struct): Options from the
%         input es struct are used with a new initial value for the
%         pseudo random number generator. Calling
%         es_options_new('es', es) is identical with calling
%         es_options_new('options', es_get('options', es)). 
%
% The evolution strategy is controlled through the options
% struct. In particular the random seed is set in
% es_options_new. To replicate results use exactly the same options
% struct. To modify an existing options struct use es_options_set.
%
% Examples: options=es_options_new;   
%                or changing a default value with
%           options=es_options_new('TolFun', 1e-5); 
%                which is identical with
%           options=es_options_set(es_options_new, 'TolFun', 1e-5); 
%                Or calling
%           opts=es_options_new('MaxFunEvals', 4000);
%                which yields in case of problem dimension 20
%                the same result as
%           opts=es_options_new('MaxFunEvals', '10*N^2');
%                Or with a longer list
%           opts=es_options_new('TolFun', 1e-5, 'TolX', 1e-4); 
%                To set these options in a new es struct do
%           es = es_new(N, xlow, xup, opts);
%                and in an existing one do
%           es = es_set(es, 'options', opts);
%
% For further examples call es_help (or see file es_READ.ME). 
%
% See also es_options_set, es_new, es_get, es_set. 

%--------%--------%--------%--------%--------%--------%--------%--------%
  % Input arguments will be treated calling es_options_set at the end

  % Stop criteria
  options.maxfunevals = inf;
  options.tolfun = 0;
  options.tolx = 1e-7; % can be N by 1 matrix, stop if *all*...
  options.minfunvalue = -inf;
  options.numprecstop = 0; % flag: stop in case of numeric precision problem
  options.maxcondstop = 0; % flag: stop in case of maximal condition reached

  % Initial Values
  options.xstart = []; % ==> ones(N,1).*(xup+xlow)/2;
  options.sigmafacstart = 1;
  options.randomseed = 100*sum(clock);

  % Initial scaling
  options.scaling = []; % ==> ones(N,1).*(xup-xlow)/2
  
  % Bounds
  options.lbound = -inf;
  options.ubound = inf;

  % Minimal standard deviation in each coordinate axis
  options.minxchange = 0; 
			      
  % Static strategy parameters which are N-dependent
  options.lambda = '4+floor(3*log(N))';
  options.mu = 'floor(lambda/2)'; % may be set also in es_options_set
  options.recoweights = 'log((lambda+1)/2)-log(1:mu)'; % s.a.
  options.cc = '4/(N+4)';
  options.ccov = '2/(N+2^0.5)^2';
  options.cs = '4/(N+4)';
  options.damp = '(1-min(0.7,N/maxiter))/cs + 1'; 

  % Miscellaneous 
  options.initphase = 1; % flag: wait with adaptation of the
                         % covariance matrix until step-size has
                         % become sufficient large 
  options.updatemod = 1; % update C if mod(countgen, updatemod) < 1
  
  if ~isempty(varargin)
    options = es_options_set(options, varargin{:});
  end
%--------%--------%--------% END es_options_new %--------%--------%--------%
