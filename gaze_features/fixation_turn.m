function [fix_rate, mean_fix_dur] = fixation_turn(fix_start_time, fix_end_time, turn_start, turn_end)
% Computes fixation statistics per turn interval
%
% INPUTS:
%   fix_start_time : nFix × 1 vector of fixation onset times (seconds)
%   fix_end_time   : nFix × 1 vector of fixation offset times (seconds)
%   turn_start     : nTurn × 1 vector of turn start times (seconds)
%   turn_end       : nTurn × 1 vector of turn end times (seconds)
%
% OUTPUTS:
%   fix_rate       : nTurn × 1 vector
%                    fixations per second during each turn
%
%   mean_fix_dur   : nTurn × 1 vector
%                    mean fixation duration (seconds) during each turn
%                    NaN if no fixations occur in the turn
%
% NOTES:
%   - A fixation is counted if it overlaps the turn interval
arguments
    fix_start_time double
    fix_end_time double
    turn_start double
    turn_end double
end

%% Init
fix_start_time = fix_start_time(:);
fix_end_time   = fix_end_time(:);
turn_start     = turn_start(:);
turn_end       = turn_end(:);

nTurn = length(turn_start);

% Sanity checks
if length(fix_start_time) ~= length(fix_end_time)
    error('Fixation start/end vectors must have the same length');
end

if length(turn_start) ~= length(turn_end)
    error('turn_start and turn_end must have the same length');
end

fix_rate     = zeros(nTurn,1);
mean_fix_dur = NaN(nTurn,1);

%% Loop over turns
for t = 1:nTurn
    
    % Duration of the current turn
    turn_duration = turn_end(t) - turn_start(t);
    
    if turn_duration <= 0
        fix_rate(t) = NaN;
        continue
    end
    
    % Find fixations overlapping this turn
    idx = fix_start_time < turn_end(t) & ...
          fix_end_time   > turn_start(t);
        
    % Fixation rate (fixations per second)
    fix_rate(t) = sum(idx) / turn_duration;
    
    % Mean fixation duration
    if n_fix > 0
        durations = fix_end_time(idx) - fix_start_time(idx);
        mean_fix_dur(t) = mean(durations, 'omitnan');
    end
end

end
