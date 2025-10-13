function [VAD, vThreshold] = levels2vad( vL, fs, varargin )
% levels2vad - Compute Voice Activity Detection (VAD) from level data
%
% Syntax:
%   [VAD, vThreshold] = levels2vad(vL, fs, varargin)
%
% Description:
%   This function computes a Voice Activity Detection (VAD) signal from
%   input level data. It uses clustering and smoothing to determine
%   which parts of the signal are likely to be speech.
%
% Inputs:
%   vL          - Input level data. Can be a numeric array or a string
%                 (see Note).
%   fs          - Sampling rate of levels in Hz
%   varargin    - Key-value pairs for configuration parameters; see
%                 "levels2vad help" for details
%
% Outputs:
%   VAD         - Binary voice activity detection signal (logical array).
%
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
  sCfg.smooth = 0.09;
  sHelp.smooth = 'smoothing time of decision output, for short gap/short speech removal';
  sCfg.hannwnd = 1;
  sHelp.hannwnd = 'duration of von-Hann window for filtering in intensity domain';
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
  len_hannwnd = round(fs*sCfg.hannwnd);
  len_smooth = round(fs*sCfg.smooth);
  vThreshold = zeros(1,num_channels);
  for ch=1:num_channels
    levels = vL(:,ch);
    
    if len_hannwnd > 0
      levels = smoothed_levels( levels, len_hannwnd );
    end
    levels(levels<sCfg.minlevel) = [];
    l_range = quantile(levels,1-[sCfg.maxcontrib,sCfg.mincontrib]);
    idx = kmeans( levels, 2 );
    v_means = sort([mean(levels(idx==1)),mean(levels(idx==2))]);
    threshold = mean(v_means);
    threshold = max(threshold, l_range(1));
    threshold = min(threshold, l_range(2));
    VAD(:,ch) = vL(:,ch) > threshold;
    vThreshold(ch) = threshold;
  end
  if len_smooth > 0
    VAD = filtfilt( ones(len_smooth,1)/len_smooth, 1, VAD ) >= 0.5;
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

function l = smoothed_levels( l, winlen )
% smoothed_levels - calculated smoothed levels from short-term
%                   logarithmic levels
%
% l = smoothed_levels( l, winlen )
%
% l      : short-term logarithmic (dB) levels
% winlen : window length in blocks
%
% The return value is the level, smoothed with a von-Hann window of winlen
% length. The lsmooth is time-aligned to correspond to a von-Hann
% window symmetric around zero.

% test for correct input data:
    if (numel(size(l)) ~= 2)||(min(size(l)) ~= 1)||(prod(size(l))==1)
        error('l must be a vector');
    end
    if numel(l) <= winlen
        error('winlen must be smaller than length of l');
    end
    % generate von-Hann window:
    win = hann( winlen );
    % normalize window for enegry preservation:
    win = win / sum(win);
    % convert dB levels to MS levels:
    l = 10.^(0.1*l(:));
    padlen = floor(winlen/2);
    % apply filter
    lsmooth = fftfilt(win,[zeros(padlen,1)+l(1);l;zeros(winlen,1)+l(end)]);
    % select relevant portion:
    l = lsmooth(winlen-1+(1:numel(l)));
    % convert back to dB levels:
    l = 10*log10(l);
end
