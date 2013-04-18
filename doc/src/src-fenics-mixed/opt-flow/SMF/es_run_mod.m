function es = es_run_mod(es, strfun, niter, param)
% function es = es_run(es, strfun [, niter])
%   minimizes the function named strfun for maximal niter iterations
%
% Input Arguments:
%   es : evolution strategy struct, see function es_new
%
%   strfun : Name of function to be minimized.  feval(strfun, X)
%            must return a scalar, where X is a N by
%            1 sized matrix, where N is the problem dimension (specified
%            with es_new(N, ...)).
%
%   niter : (optional) Number of iterations to be done at
%            most. Default is inf. The minimum number of iterations
%            is one, independ of niter. Note that niter is not
%            identical with the number of function
%            evaluations. Confer to options field lambda to find
%            out the number function evaluations per iteration. 
%
% For examples call es_help (or see file es_READ.ME). 
%
% See also es_new, es_get, es_set, es_options_new. 

% Notation: ar... is a row-array which can contain vectors or
% scalars. 

% Die Schnittstelle nach aussen sind auf den Rand korrigierte
% Punkte. 

% Test: 
%   o random fct OK
%   o Kugel OK
%   o cigar OK
%   o tablet OK
%   o elli OK
%   o twoaxes OK
%   o rosen OK: schwierig wegen der Korrelationen
%   o parab100 OK
%   o sharpR OK
%   o diffpow OK
%--------%--------%--------% BEGIN es_run %--------%--------%--------%--------%
  %% Check errors
  %%error(nargchk(2,3,nargin));
  %%if nargin < 3 
  %%  niter = inf;
  %%end
  N = es.input.N;

  %%% Display information once
  if es.countgen == 0
    % Generate string for display
    str = [];
    if es.sp.mu > 1
      if es.sp.recoweights(1) == es.sp.recoweights(end) 
	str = '_I'; 
      else 
	str = '_W';
      end
    end
    disp(['(' num2str(es.sp.mu) str ',' num2str(es.sp.lambda) ...
	  ')-CMA-ES (w=[' num2str(es.const.recow','%5.2f') '])' ]);
  end % if es.countgen == 0
  
  %%% Generation loop
  es.arflgstop = [];
  stopgen = es.countgen + niter;
  while isempty(es.arflgstop) & es.countgen < stopgen
    es.countgen = es.countgen + 1;

    % Generate offspring ---modify for es_conrun---
    arz = randn(N, es.sp.lambda);
    arx = es.xreco*ones(1, es.sp.lambda) + es.sigma * (es.BD * arz);
      
    for k=1:es.sp.lambda
      % Generate valid point from offspring 
      xvalid = x_makevalid(arx(:,k), es.opts);
      
      % Calculate difference to valid point
      ardifftobound(:,k) = xvalid-arx(:,k);
      
      % Evaluate fitness of valid point
      es.arfunval(k) = feval(strfun, xvalid, param);
      es.counteval = es.counteval + 1;
    end
    
    arfunvalsorted = sort(es.arfunval);
    
    % Calculate/modify weights for bounds penalty 
    es.stbndweig = bndweig_adjust(es.stbndweig, ...
	sum(abs(ardifftobound)==0, 2)/es.sp.lambda, ...
	arfunvalsorted, es.sigma^2*diag(es.C));

    % Calculate function value with penalty term 
    es.arfunvalpen = es.arfunval + ...
	es.stbndweig.w' * ardifftobound.^2;
    
    % ---calculate and add contraints penalty for es_conrun---

    % Selection ---modify for es_conrun---
    [arfunvalpensorted, idx] = sort(es.arfunvalpen); % minimization
    xold = es.xreco;
    es.xreco = arx(:,idx(1:es.sp.mu))*es.const.recow;
    zreco = arz(:,idx(1:es.sp.mu))*es.const.recow;

    % Adapt covariance matrix ---modify for es_conrun---
    es.pc = (1-es.sp.cc)*es.pc + ...
	(es.const.ccu*es.const.cw/es.sigma) * (es.xreco-xold);
    if ~es.flginitphase
      es.C = (1-es.sp.ccov)*es.C + es.sp.ccov*es.pc*transpose(es.pc);
    end
    % adapt sigma
    es.ps = (1-es.sp.cs)*es.ps + (es.const.csu*es.const.cw) * (es.B*zreco);
    es.sigma = es.sigma * ...
	exp((norm(es.ps)-es.const.chiN)/es.const.chiN/es.sp.damp);
    % es.sigma = es.sigma * ...
    %	* exp((norm(es.ps)^2-N)/N/2/es.sp.damp);
    
    % Adjust minimal step size (es.C and es.sigma may be changed)
    idx = find(es.sigma*sqrt(diag(es.C)) < es.opts.minxchange);
    if ~isempty(idx)
      if 0
	es.sigma = 1.4*es.sigma; 
      elseif 0
	for j=idx
	  es.C(j,j) = es.C(j,j)*1.4;
	end
      elseif 1
	for j=idx
	  es.C(j,j) = (1+es.sp.ccov) * es.C(j,j);
	end
	es.sigma = max(es.sigma, max(es.opts.minxchange ./ sqrt(diag(es.C))));
      end
      % warning(['minxchange reached in coordinate(s) ' num2str(idx')]);
    end
    
    % Update B and D from C
    if mod(es.countgen, es.opts.updatemod) < 1
      es.C=triu(es.C)+transpose(triu(es.C,1)); % enforce symmetry
      [es.B,es.D] = eig(es.C);
      % limit condition of C to 1e14 + 1
      if max(diag(es.D)) > 1e14*min(diag(es.D)) 
	if es.opts.flgmaxcond > 0
	  es.arflgstop = [es.arflgstop 6];
	  warning(['maximal condition number reached' ...
		': stopflag is set']);
	elseif es.opts.flgmaxcond == 0
	  warning(['maximal condition number reached ' ...
	      '(MaxCondStop==0)']);
	  tmp = max(diag(es.D))/1e14 - min(diag(es.D));
	  es.C = es.C + tmp*eye(N); es.D = es.D + tmp*eye(N); 
	elseif es.opts.flgmaxcond == -1
	  tmp = max(diag(es.D))/1e14 - min(diag(es.D));
	  es.C = es.C + tmp*eye(N); es.D = es.D + tmp*eye(N); 
	elseif es.opts.flgmaxcond == -2
	  warning('Condition number exceeds 1e14 (MaxCondStop==-2)');
	end
      end
      es.D = diag(sqrt(diag(es.D))); % D contains standard deviations now
      es.BD = es.B*es.D; % for speed up only
    end % if mod

    % Check for numeric precision problems ---modify for es_conrun---
    if arfunvalpensorted(1) == arfunvalpensorted(1+floor(es.sp.lambda/2)) ...
          | es.xreco == es.xreco ...
	  + 0.2*es.sigma*es.BD(:,1+floor(mod(es.countgen,N)))
      tstr = '';
      if arfunvalpensorted(1) == arfunvalpensorted(1+floor(es.sp.lambda/2))
	tstr = ' (function values equal)';
      end
      if es.opts.flgnumprec > 0
	es.arflgstop = [es.arflgstop 5];
	warning(['There seems to be a numeric precision problem' ...
	      tstr ': stopflag is set']);
      elseif es.opts.flgnumprec == 0
	es.sigma = 1.4*es.sigma; 
	warning(['There seems to be a numeric precision problem' ...
	      tstr '): sigma is increased']);
      elseif es.opts.flgnumprec == -1
	es.sigma = 1.4*es.sigma; 
      elseif es.opts.flgnumprec == -2
	warning(['There seems to be a numeric precision problem' tstr]);
      end
    end
    
    % Test for end of initial phase
    if es.flginitphase & es.countgen > 2/es.sp.cs
      if (norm(es.ps)-es.const.chiN)/es.const.chiN < 0.05 
	es.flginitphase = 0;
      end
    end
    
    % Check stop criteria
    if es.counteval >= es.opts.maxfunevals
      es.arflgstop = [es.arflgstop 1];
    end 
    if (es.sigma / es.const.cw) * es.pc < es.opts.tolx ...
	  & es.sigma*sqrt(diag(es.C)) < es.opts.tolx
      es.arflgstop = [es.arflgstop 2];
    end 
    es.ardeltafunvalhist(1+mod(es.countgen, 3)) = ...
	arfunvalsorted(end) - arfunvalsorted(1);
    if max(es.ardeltafunvalhist) < es.opts.tolfun
      es.arflgstop = [es.arflgstop 3];
    end
    if arfunvalsorted(1) <= es.opts.minfunvalue 
      es.arflgstop = [es.arflgstop 4];
    end
    
    %%% Record data
    [es.rec.funval, idx] = min(es.arfunval);
    es.rec.x = x_makevalid(arx(:,idx), es.opts);

    % save overall best point
    if es.countgen == 1 
      es.rec.best.funval = es.rec.funval;
      es.rec.best.x = es.rec.x;
      es.rec.best.eval = es.counteval;
    elseif es.rec.funval < es.rec.best.funval
      es.rec.best.funval = es.rec.funval;
      es.rec.best.x = es.rec.x;
      es.rec.best.eval = es.counteval;
    end
    
  end % while isempty(es.arflgstop)
%--------%--------%--------% END es_run %--------%--------%--------%--------%

%***************************************************************************%
function st = bndweig_adjust(st, feasratio, arfunval, sigsqrdiagC)
% function stw = bndweig_adjust(stw, feasratio, arfunval, sigsqrdiagC)
%   modifies the weights for bounds in stw
%
% Input:
%   stw : weights-struct which may be empty ("first call") or is
%     the output of former calls of bndweig_adjust. stw contains
%     cumrf1 : desired minimal cumulative ratio of feasible
%       solutions in time1 generations.
%     cumrf2 : desired maximal cumulative ratio of feasible
%       solutions in time2 generations.
%     time1 : Number of generations for cumrf1.
%     time2 : Number of generations for cumrf2.
%     count : Number of former calls.
%     w : weights, Nx1-vector.
%     lenrfrec : Length of array of recorded ratios of feasible
%     arrfrec : Array of ratios of feasible, recorded for lenrfrec
%       generations.
%     timelastset : generation count, when w were set last time.
%     flglargercumrf2 : holds result from sum(...) > cumrf2 of
%                       last generation
%     ardeltq : function value differences of a few generations. 
%     ideltq : Index for ardeltq.
%
%   feasratio   : Nx1-vector of the ratio of feasible points in the
%                 current population/generation
%   arfunval    : sorted 1xlambda-array of function values
%   sigsqrdiagC : es.sigma^2*diag(es.C)
%
%

%   Ziel: Im Fall einer "aktiven" Grenze soll fuer den Anteil der gueltigen
%   Punkte rf in jeder Koordinate, in jeder Generation, gelten:
%
%         cumrf1/time1 < rf < cumrf2/time2
%
%     mit
%       cumrf1 = 0.1, time1 = 3*N+9, cumrf2 = 0.5, time2 = 1
%
%   Algorithmus: 
%   1) Falls die Summe ueber time2 den Wert cumrf2 "kreuzt", d.h.
%         sum_{t = g-time2  ... g } rf^(t) > cumrf2   und
%         sum_{t = g-time2-1...g-1} rf^(t) < cumrf2, oder umgekehrt, 
%      und das letzte Setzen der Gewichte mindestens N/4+4
%      Generationen her ist, setze
%         w = deltaQ/N/sigma^2/diag(C)
%   2) Falls sum_{t=g-time1...g} rf^(t) < cumrf1
%      vergroessere w. 
%
%   Diskussion: Fuer rf > cumrf2/time2 liegt das Minimum vermutlich in
%   der gueltigen Region. Wenn nicht, ist die Formel fuer das Setzen
%   der Gewichte "buggy", d.h. es wurde eine Situation nicht
%   beruecksichtigt. Gegen eine Verkleinerung von Gewichten spricht
%   a) die Gewicht koennten sehr klein werden und die Strategie
%   weit in den ungueltigen Bereich laufen und b) die Adaptation
%   der Mutationsverteilung wird sich stetig an die Verkleinerung
%   adaptieren. Das erschwert die Interpretation des
%   Strategieverhaltens. 
%
%--------%--------%--------% BEGIN bndweig_adjust %--------%--------%--------%
  N = size(feasratio, 1);
  % instead of function bndweig_new
  if isempty(st) 
    st.count = 0;
    st.cumrf1 = 0.1;
    st.time1 = 3*N+9;
    st.cumrf2 = 0.5;
    st.time2 = 1;
    st.w = zeros(N,1);
    st.lenrfrec = max(st.time1, st.time2);
    st.arrfrec = ones(N, st.lenrfrec);
    st.timelastset = -inf*ones(N,1); 
    st.flglargercumrf2 = ones(N,1);
    st.ardeltq = [];
    st.ideltq = 0;
    st.out1 = []; st.out2 = []; st.out3 = []; st.out4 = [];
    st.outt0=clock; st.outetime=0.0; 
  elseif N ~= size(st.w, 1)
    error('bug detected');
  end

  st.count = st.count+1;

  %%% Manage ardeltq
  lambda = size(arfunval, 2);
  dq = arfunval(ceil(3*lambda/4)) - arfunval(ceil(lambda/4));
  % hold dq for a few generations
  if dq > 0
    % ++st.ideltq :
    st.ideltq = mod(st.ideltq, 10) + 1;
    st.ardeltq(st.ideltq) = dq;
  end
  
  if isempty(st.ardeltq)
    error(['Function values were constant. Starting point was probably' ...
	  ' beyond the given bounds or initial search cubus widely exceeds' ...
	  ' the given bounds.']);
  end
  
  %%% Update arrfrec
  irfrec = mod(st.count-1, st.lenrfrec) + 1;
  st.arrfrec(:, irfrec) = feasratio;
  
  flgTEST = 0;
  
  %%% Operations for count <= 3
  if st.count < 3
    return;
  elseif flgTEST & st.count == 3
    st.w = (min(st.ardeltq)/N/3) ./ sigsqrdiagC;
    st.timelastset(:) = st.count;
    return;
  end
  
  %%% Enlarge weights
  % set idx for summation over time1 data
  idx = mod((irfrec-st.time1):(irfrec-1), st.lenrfrec) + 1;
  % find components, where sum(rf) < cumrf1/time1
  idx = find(sum(st.arrfrec(:, idx), 2) < st.cumrf1);
  st.w(idx) = st.w(idx) * exp(1/(3*N+1));
  
  %%% Set weights, if sum(rf) crosses cumrf2/time2
  % set idx for summation over time2 data
  idx = mod((irfrec-st.time2):(irfrec-1), st.lenrfrec) + 1;
  % find components, where sum(rf) > cumrf2/time2
  flg = (sum(st.arrfrec(:, idx), 2) > st.cumrf2);
  % XOR of flg and st.flglargercumrf2
  idx = find(flg + st.flglargercumrf2 == 1);
  st.flglargercumrf2 = flg;
  % Set not before N/4+4 generations past
  idx = intersect(idx, find(st.count-st.timelastset > N/4+4));
  st.w(idx) = (min(st.ardeltq)/N/3) ./ sigsqrdiagC(idx);
  st.timelastset(idx) = st.count;
  
  if flgTEST
  st.out1 = [st.out1; feasratio'];
  st.out2 = [st.out2; ...
	(sum(st.arrfrec(:,mod((st.count-ceil(N+2)):(st.count-1), ...
	st.lenrfrec)+1), 2)/ceil(N+2))'];
  st.out3 = [st.out3; log10(st.w+1e-3)'];
  if length(st.ardeltq) < 10
    st.out4 = [st.out4; ones(1,10)];
  else
    st.out4 = [st.out4; log10(sort(st.ardeltq))];
  end
  if etime(clock,st.outt0) > 10*st.outetime | ...
	etime(clock, st.outt0) <= 0 % bug fix for etime
    st.outt0 = clock;
    figure(531);
    subplot(2,2,1);plot(1:size(st.out1,1), st.out1);
    title('Ratio of Feasible Solutions');grid on;zoom on; 
    subplot(2,2,2);plot(1:size(st.out2,1), st.out2);
    title('Mooving average of feasible ratios');grid on;zoom on; 
    subplot(2,2,3);plot(1:size(st.out3,1), st.out3);
    title('Weights');grid on;zoom on; 
    subplot(2,2,4);plot(1:size(st.out4,1), st.out4);
    title('deltaQ, sorted');grid on;zoom on; 
    st.outetime = etime(clock, st.outt0);
  end
  end
%--------%--------%--------% END bndweig_adjust %--------%--------%--------%

  
function x=x_makevalid( x, opts)
  idx = find(x < opts.lbound);
  x(idx) = opts.lbound(idx);
  idx = find(x > opts.ubound);
  x(idx) = opts.ubound(idx);
  
