function VAD = levels2vad( vL, varargin )
% levels2vad - Compute Voice Activity Detection (VAD) from level data
%
% Syntax:
%   [VAD, VAD2] = levels2vad(vL, varargin)
%
% Description:
%   This function computes a Voice Activity Detection (VAD) signal from
%   input level data. It uses clustering and smoothing to determine
%   which parts of the signal are likely to be speech.
%
% Inputs:
%   vL          - Input level data. Can be a numeric array or a string
%                 (see Note).
%   varargin    - Key-value pairs for configuration parameters; see
%                 "levels2vad help" for details
%
% Outputs:
%   VAD         - Binary voice activity detection signal (logical array).
%
% Configuration Parameters:
%   The following parameters can be specified as key-value pairs in
%   varargin:
%   - 'minlevel'  - Minimum level in dB SPL. Levels below this value
%                  are ignored in threshold estimation. Default: 25.
%   - 'mincontrib' - Minimum speech contribution. Default: 0.1.
%   - 'maxcontrib' - Maximum speech contribution. Default: 0.8.
%   - 'smooth'    - Number of level samples to smooth over for noise
%                  reduction. Default: 9.
%
% Notes:
%   - If vL is a string, it is treated as the first key in varargin.
%   - The function automatically loads required Octave packages if
%     necessary.
%   - The function uses k-means clustering to determine speech
%     thresholds.
%
% Example:
%   % Compute VAD with default parameters
%   VAD = levels2vad(levelData);
%
%   % Compute VAD with custom parameters
%   VAD = levels2vad(levelData, 'minlevel', 30, 'smooth', 5);
%
  if ischar(vL)
    varargin = {vL,varargin{:}};
  end
  sCfg = struct;
  sHelp = struct;
  sCfg.minlevel = 25;
  sHelp.minlevel = 'minimum level in dB SPL to; levels below will be ignored in threshold estimation';
  sCfg.mincontrib = 0.1;
  sHelp.mincontrib = 'minimum speech contribution';
  sCfg.maxcontrib = 0.8;
  sHelp.maxcontrib = 'maximum speech contribution';
  sCfg.smooth = 9;
  sHelp.smooth = 'number of level samples to smooth over, for gap/short speech removal';
  sCfg = parse_keyval( sCfg, sHelp, varargin{:} );
  if isempty(sCfg)
    return;
  end

  if isoctave()
    if isempty(which('kmeans'))
      pkg load statistics
    end
    if isempty(which('filtfilt'))
      pkg load signal
    end
  end
  num_channels = size(vL,2);
  VAD = zeros(size(vL));
  for ch=1:num_channels
    levels = vL(:,ch);
    levels(levels<sCfg.minlevel) = [];
    l_range = quantile(levels,1-[sCfg.maxcontrib,sCfg.mincontrib]);
    idx = kmeans( levels, 2 );
    v_means = sort([mean(levels(idx==1)),mean(levels(idx==2))]);
    threshold = mean(v_means);
    threshold = max(threshold, l_range(1));
    threshold = min(threshold, l_range(2));
    VAD(:,ch) = vL(:,ch) > threshold;
  end
  if sCfg.smooth > 0
    VAD = filtfilt( ones(sCfg.smooth,1)/sCfg.smooth, 1, VAD ) >= 0.5;
  end
  VAD = logical(VAD);
end

function b = isoctave()
  b = ~isempty(ver('Octave'));
end

function sCfg = parse_keyval( sCfg, sHelp, varargin )
  if (numel(varargin) == 1) && strcmp(varargin{1},'help')
    sOut = sprintf(' List of valid keys:\n');
    fields = fieldnames(sCfg);
    for k=1:numel(fields)
      field = fields{k};
      helps = '';
      if isfield(sHelp,field)
	helps = sHelp.(field);
      end
      defval = '?';
      if isnumeric(sCfg.(field)) || islogical(sCfg.(field))
	defval = mat2str(sCfg.(field));
      end
      if ischar(sCfg.(field))
	defval = sCfg.(field);
      end
      if isa(sCfg.(field),'function_handle')
	defval = ['@',func2str(sCfg.(field)),' (function handle)'];
      end
      if iscellstr(sCfg.(field))
	defval = '{';
	if isempty(sCfg.(field))
	  defval(end+1) = ' ';
	end
	for ke=1:numel(sCfg.(field))
	  defval = sprintf('%s''%s'',',defval,sCfg.(field){ke});
	end
	defval(end) = '}';
      end
      sOut = sprintf('%s - %s:\n   %s\n   (default: %s)\n\n',...
		     sOut,field, ...
		     helps,defval);
    end
    disp(sOut);
    sCfg = [];
    return
  end
  if isempty(varargin)
    return;
  end
  if mod(numel(varargin),2) > 0
    error('Invalid (odd) number of input arguments.');
  end
  for k=1:2:numel(varargin)
    if ~ischar(varargin{k})
      error(sprintf('key %d is not a string.',k));
    end
    field = lower(varargin{k});
    if ~isfield(sCfg,field)
      error(sprintf('key %d (''%s'') is not a valid key.',k, ...
		    field));
    end
    sCfg.(field) = varargin{k+1};
  end
end