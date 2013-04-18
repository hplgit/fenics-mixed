function opts=es_options_set(opts, varargin)
% function opts=es_options_set(opts[, 'param1', value1, 'param2', value2,...])
%   sets new values in the options struct opts
%
% Input arguments: 
%   opts : Options struct, see es_options_new
%   varargin : See es_options_new.
%
% For examples call es_help (or see file es_READ.ME). 
%
% See also es_options_new, es_new

% Function es_options_set only sets fields in the opts
% struct. These fields are evaluated in function es_set. 
%--------%--------%--------% BEGIN es_options_set %--------%--------%--------%
  flgmu = 0; flgrecow = 0; args = varargin;
  for i = 1:2:length(args)
    if length(args) < i+1
      error(['Value for ' args{i} ' missing']);
    end
    flgdone = 0;
    if i == 1 % check for struct es or options 
      if strcmp('es', args{1})
	opts = es_options_set(opts, 'options', es_get(args{2}, 'options'));
	flgdone = 1;
      elseif strcmp('options', args{1})
	opts = args{2};
	opts.randomseed = 100*sum(clock);
	flgdone = 1;
      end
    end
    if ~flgdone % not es or options
      if ~isfield(opts, lower(args{i}))
	error([lower(args{i}) ' is not a field']);
      end
      opts = setfield(opts, lower(args{i}), args{i+1});
      if strcmp('mu', lower(args{i}))
	if ~flgrecow
	  opts.recoweights = 'ones(1,mu)';
	end
      elseif strcmp('recoweights', lower(args{i}))
	if isempty(opts.recoweights)
	  error('recoweights cannot be set to []');
	end
	flgrecow = 1;
      end
    end
  end % for
  
%--------%--------%--------% END es_options_set %--------%--------%--------%
