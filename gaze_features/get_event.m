function onset_event = get_event(gaze_in_roi, time, timestamps)
% Returns the ROI or event(e.g. blink, fixation,...) active at given timestamps
%
% INPUTS:
%   gaze_in_roi : nTime × nROI matrix OR nTime × 1 vector
%                 - ROI case: gaze_in_roi(t,r) = 1 if gaze is in ROI r at time t
%                 - Event case (e.g. blink): gaze_in_roi(t) = 1 if event active
%
%   time        : nTime × 1 vector of timestamps (seconds)
%   timestamps  : nEvent × 1 vector of event times (e.g. turn_start / turn_end)
%
% OUTPUT:
%   onset_event : nEvent × 1 vector
%                 - ROI index (1..nROI) at each timestamp
%                 - 1 for binary event present
%                 - 0 if no ROI / event present
%
% NOTES:
%   - Nearest-neighbor matching is used between time and timestamps
%   - If multiple ROIs are active, the first ROI index is returned
arguments
    gaze_in_roi logical
    time double
    timestamps double
end


%% Init
time = time(:);
timestamps = timestamps(:);

if size(gaze_in_roi,1) ~= length(time)
    error('gaze_in_roi and time must have the same number of rows');
end

nTimestamps = length(timestamps);
nROIs = size(gaze_in_roi,2);

% Ensure logical
gaze_in_roi = logical(gaze_in_roi);


onset_event = zeros(nTimestamps,1);

%%  Loop over timestamps
for i = 1:nTimestamps
    
    % Find nearest time index
    [~, idx] = min(abs(time - timestamps(i)));
    
    % Check event / ROI state at that time
    if nROIs == 1
        % Binary event (e.g. blink)
        onset_event(i) = gaze_in_roi(idx);
    else
        % ROI case
        roi_idx = find(gaze_in_roi(idx,:), 1, 'first');
        if ~isempty(roi_idx)
            onset_event(i) = roi_idx;
        end
    end
end

end
