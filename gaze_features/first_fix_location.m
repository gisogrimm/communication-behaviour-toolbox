function first_fix_loc = first_fix_location(fix_in_rois, fix_start_time, fix_end_time, turn_start, turn_end)
% Determines the ROI of the first fixation after a turn
%
% INPUTS:
%   fix_in_rois    : nFix × nROIs logical or numeric matrix
%                    fix_in_rois(i,j) = 1 if fixation i is in ROI j
%
%   fix_start_time : nFix × 1 vector of fixation start times (seconds)
%   fix_end_time   : nFix × 1 vector of fixation end times (seconds)
%   turn_start     : nTurn × 1 vector of turn onset times (seconds)
%   turn_end     : nTurn × 1 vector of turn offset times (seconds)

%
% OUTPUT:
%   first_fix_loc  : nTurn × 1 vector
%                    ROI index (1..nROIs) of first fixation after each turn
%                    NaN if no fixation occurs after the turn
%
% NOTES:
%   - A fixation is considered "after a turn" if its start time is >= turn_start
%   - If a fixation belongs to multiple ROIs, the first matching ROI index
%     is returned
%
arguments
    fix_in_rois logical
    fix_start_time double
    fix_end_time double
    turn_start double
    turn_end double
end

%% Init
fix_start_time = fix_start_time(:);
fix_end_time   = fix_end_time(:);
turn_start     = turn_start(:);
turn_end     = turn_end(:);

[nFix, ~] = size(fix_in_rois);
nTurn = length(turn_start);

% Sanity checks
if length(fix_start_time) ~= nFix || length(fix_end_time) ~= nFix
    error('fix_in_rois and fixation time vectors must have compatible sizes');
end

if length(turn_end) ~= nTurn
    error('turn_start and turn_end vectors must have compatible sizes');
end

first_fix_loc = NaN(nTurn, 1);

%% Find fix
% Loop over turns
for t = 1:nTurn
    
    % Find fixations starting after the turn
    idx_fix = find(fix_end_time >= turn_start(t), 1, 'first');
    
    % If no fixation found or if fixation starts after turn end, leave NaN
    if isempty(idx_fix)
        continue
    end
    
    if fix_start_time(idx_fix)>turn_end(t)
        continue
    end
    
    % Find ROI(s) of this fixation
    roi_idx = find(fix_in_rois(idx_fix, :), 1, 'first');
    
    % Assign ROI index if fixation is in any ROI
    if ~isempty(roi_idx)
        first_fix_loc(t) = roi_idx;
    end
end

end
