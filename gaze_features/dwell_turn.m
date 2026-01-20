function [mean_dwell_time, total_dwell_time, revisit_rate] = ...
    dwell_turn(dwell_start, dwell_end, turn_start, turn_end)
% Computes dwell statistics per turn for a single ROI
%
% INPUTS:
%   dwell_start : nDwell × 1 vector of dwell onset times (seconds)
%   dwell_end   : nDwell × 1 vector of dwell offset times (seconds)
%   turn_start  : nTurn × 1 vector of turn start times (seconds)
%   turn_end    : nTurn × 1 vector of turn end times (seconds)
%
% OUTPUTS:
%   mean_dwell_time  : nTurn × 1 vector
%                      mean dwell duration (seconds) during each turn
%                      NaN if no dwells occur in the turn
%
%   total_dwell_time : nTurn × 1 vector
%                      total dwell time (seconds) during each turn
%
%   revisit_rate     : nTurn × 1 vector
%                      number of dwell onsets per second during each turn
%
% NOTES:
%   - A dwell is counted if it overlaps the turn interval
arguments
    dwell_start double
    dwell_end double
    turn_start double
    turn_end double
end

%% Init

dwell_start = dwell_start(:);
dwell_end   = dwell_end(:);
turn_start  = turn_start(:);
turn_end    = turn_end(:);

nTurn = length(turn_start);

% Sanity checks
if length(dwell_start) ~= length(dwell_end)
    error('dwell_start and dwell_end must have the same length');
end

if length(turn_start) ~= length(turn_end)
    error('turn_start and turn_end must have the same length');
end

mean_dwell_time  = NaN(nTurn,1);
total_dwell_time = zeros(nTurn,1);
revisit_rate     = zeros(nTurn,1);

%% calc dwell features
%Loop over turns
for t = 1:nTurn
    
    % Duration of this turn
    turn_duration = turn_end(t) - turn_start(t);
    
    if turn_duration <= 0
        mean_dwell_time(t)  = NaN;
        total_dwell_time(t) = NaN;
        revisit_rate(t)     = NaN;
        continue
    end
    
    % Find dwells overlapping this turn
    idx = dwell_start < turn_end(t) & dwell_end   > turn_start(t);
    
    n_dwell = sum(idx);
    
    if n_dwell == 0
        mean_dwell_time(t)  = NaN;
        total_dwell_time(t) = 0;
        revisit_rate(t)     = 0;
        continue
    end
    
    % Compute overlap-corrected dwell durations
    overlap_start = max(dwell_start(idx), turn_start(t));
    overlap_end   = min(dwell_end(idx),   turn_end(t));
    
    dwell_durations = overlap_end - overlap_start;
    
    % Total dwell time
    total_dwell_time(t) = sum(dwell_durations);
    
    % Mean dwell duration
    mean_dwell_time(t) = mean(dwell_durations, 'omitnan');
    
    % Revisit rate (dwells per second)
    revisit_rate(t) = n_dwell / turn_duration;
end

end
