function [fix_idx,fix_start_idx,fix_end_idx,durations,position,precision_rms,gaze_velocity] = detect_fixations(gaze,time,options)
% DETECT_FIXATIONS identifies fixation periods based on a velocity threshold.
%
% INPUTS:
%   gaze   : gaze position
%   time   : either timestamps in s  OR sampling frequency
%   options:
%       velocityThreshold : velocity threshold
%       minDuration       : minimum fixation duration (ms)
%       mergeTime         : temporal gap threshold for merging fixations (ms)
%       mergeAngle        : spatial threshold for merging fixations
%
% OUTPUTS:
%   fix_idx       : fixation index per sample (0 = no fixation)
%   fix_start_idx : starting sample index of each fixation
%   fix_end_idx   : ending sample index of each fixation
%   durations     : fixation durations (s)
%   position      : mean gaze position during each fixation
%   precision_rms : RMS noise within each fixation (spatial stability)
%

arguments
    gaze double
    time double
    options.velocityThreshold double = 100
    options.minDuration double = 0.100       % s
    options.mergeTime double = 0.050          % s
    options.mergeAngle double = 0.1        % degrees
end


%% -------------------------------------------
%  Normalize gaze and time dimensionality
%  Make sure gaze is n×d (not d×n)
% --------------------------------------------
[r, c] = size(gaze);

% If gaze is 1×n, 2×n, or 3×n → transpose to n×d
if r <= 3 && c >= 1
    gaze = gaze';
end

if size(time,1)>1
    time = time';
end

%% -------------------------------------------
%  Compute gaze velocity magnitude
% --------------------------------------------

% Euclidean distance between consecutive gaze samples (1D/2D/3D general)
gazeDiff = sqrt(sum(diff(gaze,1,1).^2, 2));   % (n-1)×1

% If `time` is a scalar → it's sampling frequency
if length(time) == 1
    gaze_velocity = gazeDiff .* time;         % v = dx * fs
    fs = time;
else
    % Timestamps provided → compute dt
    timeDiff = diff(time);
    gaze_velocity = gazeDiff ./ timeDiff;     % v = dx / dt
    fs = 1 / median(timeDiff);                % estimated sampling frequency
end

% Pad first sample (velocity undefined there)
gaze_velocity = [0; abs(gaze_velocity)];


%% -------------------------------------------
%  Detect fixation periods (velocity below threshold)
% --------------------------------------------
tmp_fix_idx = gaze_velocity <= options.velocityThreshold;

% Extract fixation periods using helper function
[fix_start_idx, fix_end_idx] = findMovementPeriod( ...
    tmp_fix_idx, ...
    ceil(fs * options.minDuration * 1e-3), ...   % minimum samples
    1, ...                                       % allow merging
    gaze(:,1), ...                                % position (x) for merge angle
    ceil(fs * options.mergeTime * 1e-3), ...      % merge time threshold
    options.mergeAngle);                          % merge angle threshold

% Create sample-level fixation index array
fix_idx = zeros(size(gaze_velocity));

if ~isempty(fix_start_idx)
    for i = 1:length(fix_start_idx)
        fix_idx(fix_start_idx(i):fix_end_idx(i)) = i;
    end
end


%% -------------------------------------------
%  Compute fixation metrics
% --------------------------------------------

% Durations based on timestamps
durations = time(fix_end_idx) - time(fix_start_idx);

% mean position for all available dimensions (1D, 2D, or 3D)
position = arrayfun(@(a,b) mean(gaze(a:b,:), 1), ...
    fix_start_idx, fix_end_idx, 'UniformOutput', false);
position = cell2mat(position');


% RMS precision inside each fixation
precision_rms = arrayfun(@(i) ...
    sqrt(sum(gazeDiff(fix_start_idx(i):fix_end_idx(i))) / ...
    (fix_end_idx(i) - fix_start_idx(i) + 1)), ...
    1:length(fix_start_idx));


%% ------------------------------------------------------------------------
%  Helper function: find fixation start and end indices
% ------------------------------------------------------------------------
    function [idx_starts, idx_ends] = findMovementPeriod(data,min_samples,merge,gaze_position,merge_threshold,merge_angle)
        % FINDMOVEMENTPERIOD finds the start/end indices of contiguous TRUE segments
        % in a logical array, applies merging rules, and removes short segments.

        if nargin < 2, min_samples = 1; end
        if nargin < 3, merge = 0; end
        if nargin < 5, merge_threshold = 5; end
        if nargin < 6, merge_angle = 2; end

        % Logical indices of samples meeting velocity criterion
        idx = find(data);

        if isempty(idx)
            idx_starts = [];
            idx_ends = [];
            return;
        end

        % Extract start & end indices of contiguous TRUE runs
        idx_starts = [max(idx(1),2); idx(find(diff(idx) ~= 1) + 1)] - 1;
        idx_ends   = [idx(diff(idx) ~= 1); idx(end)];

        % -----------------------------------------------
        % Merge close fixation periods
        % -----------------------------------------------
        if merge
            gaps = idx_starts(2:end) - idx_ends(1:end-1);
            gap_size = abs(gaze_position(idx_starts(2:end)) - gaze_position(idx_ends(1:end-1)));

            merge_idx = find(gaps < merge_threshold & gap_size < merge_angle);

            if ~isempty(merge_idx)
                idx_starts(merge_idx+1) = [];
                idx_ends(merge_idx) = [];
            end
        end

        % -----------------------------------------------
        % Remove very short fixations
        % -----------------------------------------------
        durations_local = idx_ends - idx_starts;
        remove_idx = find(durations_local < min_samples);

        if ~isempty(remove_idx)
            idx_starts(remove_idx) = [];
            idx_ends(remove_idx) = [];
        end
    end % helper function

end % main function
