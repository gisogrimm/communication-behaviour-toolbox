function [start_time, end_time, mean_duration, total_duration, revisits] = ...
         dwell_roi_features(in_roi, time, options)
% DWELL_ROI_FEATURES Extracts dwell/fixation features for a region of interest
%
% INPUTS:
%   in_roi     - n×1 binary vector: 1 = gaze inside ROI, 0 = outside
%   time       - n×1 vector of timestamps (seconds)
%   options.min_dwell  - minimum dwell duration in milliseconds (optional, default = 50 ms)
%   options.merge_isi  - maximum gap to merge consecutive dwells in milliseconds (optional, default = 50 ms)
%
% OUTPUTS:
%   start_time    - vector of start times for each dwell
%   end_time      - vector of end times for each dwell
%   mean_duration - mean duration of dwell events (seconds)
%   total_duration- total time spent in ROI (seconds)
%   revisits      - number of times ROI was revisited (number of dwell events - 1)

arguments
    in_roi double
    time double
    options.min_dwell double = 50    % ms
    options.merge_isi double = 50    % ms
end

%% Ensure column vectors
if size(in_roi,2) > 1, in_roi = in_roi'; end
if size(time,2) > 1, time = time'; end

%% Find start and end indices of dwells
d_in_roi = diff([0; in_roi; 0]);  % pad with zeros to catch edges
start_idx = find(d_in_roi == 1);  % rising edge
end_idx   = find(d_in_roi == -1) - 1; % falling edge

%% Convert indices to time
start_time = time(start_idx);
end_time   = time(end_idx);
durations  = end_time - start_time;


%% Merge consecutive dwells separated by less than merge_isi
i = 1;
while i < length(start_time)
    gap = start_time(i+1) - end_time(i); % seconds
    if gap*1000 <= options.merge_isi
        % Merge current and next dwell
        end_time(i) = end_time(i+1);
        start_time(i+1) = [];
        end_time(i+1) = [];
        durations(i) = end_time(i) - start_time(i);
        durations(i+1:end) = durations(i+1:end); % keep remaining durations
    else
        i = i + 1;
    end
end

%% Remove dwells shorter than min_dwell
valid = durations*1000 >= options.min_dwell; % convert s -> ms
start_time = start_time(valid);
end_time   = end_time(valid);
durations  = durations(valid);

%% Compute summary metrics
mean_duration  = mean(durations);
total_duration = sum(durations);
revisits       = max(0, length(durations)-1);

end
