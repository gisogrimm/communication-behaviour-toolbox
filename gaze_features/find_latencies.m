function latencies = find_latencies(gaze_times, turn_times, timewindow)
%computes gaze latencies relative to turn events
%
% INPUTS:
%   gaze_times  : vector of gaze event timestamps (seconds)
%   turn_times  : vector of turn event timestamps (seconds)
%   timewindow  : 
%       - scalar T     → window is [-T, +T] seconds around turn
%       - 2-element    → window is [-timewindow(1), +timewindow(2)]
%
% OUTPUT:
%   latencies   : vector of latencies (seconds), one per turn event
%                 NaN if no gaze event occurs in the time window
%
% Latency is defined as:
%   gaze_time - turn_time
%
arguments
    gaze_times double
    turn_times double
    timewindow double
end

%% Init
gaze_times = gaze_times(:);
turn_times = turn_times(:);

% Interpret time window
if isscalar(timewindow)
    t_before = timewindow;
    t_after  = timewindow;
elseif numel(timewindow) == 2
    t_before = timewindow(1);
    t_after  = timewindow(2);
else
    error('timewindow must be a scalar or a 2-element vector');
end


%% Compute latencies
latencies = NaN(size(turn_times));   % default output

for i = 1:length(turn_times)
    
    % Time window relative to current turn
    t0 = turn_times(i);
    
    % Find gaze events inside window
    idx = gaze_times >= t0 - t_before & gaze_times <= t0 + t_after;
    
    if any(idx)
        % Use nearest gaze event in time
        [~, nearest_idx] = min(abs(gaze_times(idx) - t0));
        valid_times = gaze_times(idx);
        latencies(i) = valid_times(nearest_idx) - t0;
    end
end

end
