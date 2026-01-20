function [saccade_rate, avg_peak_velocity] = saccade_turn(sacc_start_time, sacc_end_time, peak_velocity, turn_start, turn_end)
% Computes saccade statistics per turn interval
%
% INPUTS:
%   sacc_start_time : nSacc × 1 vector of saccade onset times (seconds)
%   sacc_end_time   : nSacc × 1 vector of saccade offset times (seconds)
%   peak_velocity   : nSacc × 1 vector of peak saccade velocities
%   turn_start      : nTurn × 1 vector of turn start times (seconds)
%   turn_end        : nTurn × 1 vector of turn end times (seconds)
%
% OUTPUTS:
%   saccade_rate       : nTurn × 1 vector
%                        saccades per second during each turn
%
%   avg_peak_velocity  : nTurn × 1 vector
%                        mean peak velocity of saccades during each turn
%                        NaN if no saccades occur in the turn
%
% NOTES:
%   - A saccade is counted if it overlaps the turn interval
arguments
    sacc_start_time double
    sacc_end_time double
    peak_velocity double
    turn_start double
    turn_end double
end

%% Init
sacc_start_time = sacc_start_time(:);
sacc_end_time   = sacc_end_time(:);
peak_velocity   = peak_velocity(:);
turn_start      = turn_start(:);
turn_end        = turn_end(:);

nTurn = length(turn_start);

% Sanity checks
if length(sacc_start_time) ~= length(sacc_end_time) || ...
   length(sacc_start_time) ~= length(peak_velocity)
    error('Saccade input vectors must have the same length');
end

if length(turn_start) ~= length(turn_end)
    error('turn_start and turn_end must have the same length');
end


saccade_rate      = zeros(nTurn,1);
avg_peak_velocity = NaN(nTurn,1);

%% Loop over turns
for t = 1:nTurn
    
    % Duration of this turn
    turn_duration = turn_end(t) - turn_start(t);
    
    if turn_duration <= 0
        saccade_rate(t) = NaN;
        continue
    end
    
    % Find saccades overlapping this turn
    idx = sacc_start_time < turn_end(t) & sacc_end_time   > turn_start(t);
        
    % Saccade rate (per second)
    saccade_rate(t) = sum(idx) / turn_duration;
    
    % Average peak velocity
    if n_sacc > 0
        avg_peak_velocity(t) = mean(peak_velocity(idx), 'omitnan');
    end
end

end
