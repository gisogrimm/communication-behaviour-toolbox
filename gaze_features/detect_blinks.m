function [blink_idx, blink_start_idx, blink_end_idx, durations, rate] = detect_blinks(eog, time, options)
% DETECT_BLINKS  Detects blink events in EOG data using amplitude and slope thresholds.
%
% INPUTS:
%   eog      - n×1 vertical EOG signal (microvolts).
%   time     - 1×n vector of timestamps (seconds).
%
% OPTIONS (name-value pairs):
%   options.amp_threshold      - Minimum EOG amplitude for blink detection (default: 80 µV)
%   options.slope_threshold    - Minimum derivative threshold (default: 300 µV/s)
%   options.min_blink_duration - Minimum blink duration (default: 0.05 s)
%   options.max_blink_duration - Maximum blink duration (default: 0.4 s)
%
% OUTPUTS:
%   blink_idx        - All sample indices identified as part of a blink
%   blink_start_idx  - Start index of each blink
%   blink_end_idx    - End index of each blink
%   durations        - Duration of each blink (seconds)
%   rate             - blink rate (per minute)
%
% -------------------------------------------------------------------------

arguments
    eog double
    time double
    options.amp_threshold double = 0.0008
    options.slope_threshold double = 0.003
    options.min_blink_duration double = 0.05
    options.max_blink_duration double = 0.4
end

%% Normalize dimensionality
if size(eog,1) < size(eog,2)
    eog = eog';  
end

if size(time,1) < size(time,2)
    time = time';
end

%% Compute EOG velocity (slope)
dt = diff(time);
deog = diff(eog);
velocity = deog ./ dt;

% pad to match length
velocity = [0; velocity];

%% Blink candidate detection
% Blink signature: large EOG deflection + steep slope
candidates = find(abs(eog) > options.amp_threshold | abs(velocity) > options.slope_threshold);

if isempty(candidates)
    disp("No blinks detected.");
    blink_idx = [];
    blink_start_idx = [];
    blink_end_idx = [];
    durations = [];
    rate = [];
    return
end

%% Group consecutive blink samples into events
breaks = find(diff(candidates) > 1);
blink_start_idx = [candidates(1); candidates(breaks+1)];
blink_end_idx   = [candidates(breaks); candidates(end)];

%% Compute durations
durations = time(blink_end_idx) - time(blink_start_idx);

%% Remove events based on duration
valid = durations >= options.min_blink_duration & ...
        durations <= options.max_blink_duration;

blink_start_idx = blink_start_idx(valid);
blink_end_idx   = blink_end_idx(valid);
durations       = durations(valid);
rate = 60*length(blink_start_idx)/(time(end)-time(1));

%% Return a single list of all blink sample indices
blink_idx = cell2mat(arrayfun(@(s,e) (s:e)', blink_start_idx, blink_end_idx, 'UniformOutput', false));

end
