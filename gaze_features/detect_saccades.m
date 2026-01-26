function [sacc_idx, sacc_start_idx, sacc_end_idx, durations, amplitudes, mean_velocity, peak_velocity] = detect_saccades(gaze, time, velocity, fix, options)
% DETECT_SACCADES
% Collects saccade information from gaze data using either:
%   - fixation index array (fix), or
%   - a velocity threshold (fix used as threshold)
%
% Supports 1D, 2D, or 3D gaze data.
%
% INPUTS
%   gaze      : n×d gaze samples (d = 1, 2, or 3)
%   time      : 1×n timestamps OR scalar sampling frequency
%   velocity  : gaze velocity magnitude (n×1)
%   fix       : fixation index vector OR velocity threshold
%
% OPTIONS
%   amplitude_thres     : minimum saccade amplitude (deg)
%   min_sacc_duration   : minimum saccade duration (s)
%   max_sacc_duration   : maximum saccade duration (s)
%   isi_threshold       : minimum inter-saccadic interval (ms)
%
% OUTPUTS
%   sacc_idx        : indexes that belong to a saccade
%   sacc_start_idx  : saccade onset samples
%   sacc_end_idx    : saccade offset samples
%   durations       : duration of each saccade
%   amplitudes      : Euclidean amplitude per saccade
%   mean_velocity   : mean velocity inside each saccade
%   peak_velocity   : peak velocity inside each saccade


arguments
    gaze double
    time double
    velocity double
    fix double = 50
    options.amplitude_thres double =  2
    options.min_sacc_duration double = 0.02
    options.max_sacc_duration double = 0.5
    options.isi_threshold double = 20
end


%% ----------------------------------------------------
%  Normalize gaze dimensionality (ensure n×d)
% ----------------------------------------------------
[r, c] = size(gaze);
if r <= 3 && c >= 1
    gaze = gaze';
end

% Ensure time is 1×n
if size(time,1) < 1
    time = time';
end

n = size(gaze,1);   % number of samples


%% ----------------------------------------------------
%  Determine saccade sample indices
% ----------------------------------------------------
% CASE 1: fix is an index array → define saccades as gaps in fixations
if length(fix) > 1
    % Find fixation boundaries: transitions 1→0 or 0→1
    onsets  = find(diff(fix) == -1); % fixation → non-fixation
    offsets = find(diff(fix) ==  1); % non-fixation → fixation

    % Saccades are everything not labeled as fixation
    sacc_idx = find(fix == 0);

    % Handle invalid samples (Nan gaze)
    invalid  = find(any(isnan(gaze),2));
    sacc_idx = [sacc_idx; invalid];
    
else
    % CASE 2: fix is scalar → treat as velocity threshold
    sacc_idx = find(abs(velocity) > fix);
    invalid  = find(any(isnan(gaze),2));
    sacc_idx = [sacc_idx; invalid];
end

sacc_idx = unique(sacc_idx);


%% ----------------------------------------------------
%  If no saccades found → exit
% ----------------------------------------------------
if isempty(sacc_idx)
    disp('No saccades were detected');
    sacc_start_idx=[]; sacc_end_idx=[];
    durations=[]; amplitudes=[];
    mean_velocity=[]; peak_velocity=[];
    return
end


%% ----------------------------------------------------
%  Extract saccade start and end indices (group contiguous samples)
% ----------------------------------------------------
d_idx = diff(sacc_idx);
sacc_start_idx = [sacc_idx(1); sacc_idx(find(d_idx ~= 1) + 1)];
sacc_end_idx   = [sacc_idx(find(d_idx ~= 1)); sacc_idx(end)];

% Avoid indexing past the end
sacc_end_idx(sacc_end_idx >= n) = n-1;


%% ----------------------------------------------------
%  Convert sample indices to timestamps
% ----------------------------------------------------
sacc_t_starts = time(sacc_start_idx);
sacc_t_ends   = time(sacc_end_idx);


%% ----------------------------------------------------
%  Remove saccades with too small inter-saccadic interval (ISI)
% ----------------------------------------------------
isi = sacc_t_starts(2:end) - sacc_t_ends(1:end-1);
remove_idx = find(isi < options.isi_threshold * 1e-3);

if ~isempty(remove_idx)
    sacc_start_idx(remove_idx+1) = [];
    sacc_end_idx(remove_idx)     = [];
    sacc_t_starts(remove_idx+1)  = [];
    sacc_t_ends(remove_idx)      = [];
end


%% ----------------------------------------------------
%  Compute saccade amplitudes (full 3D Euclidean)
% ----------------------------------------------------
amp_diff = gaze(sacc_end_idx,:) - gaze(sacc_start_idx,:);
amplitudes = sqrt(sum(amp_diff.^2, 2));   % n_sacc × 1


%% ----------------------------------------------------
%  Remove saccades below amplitude or duration thresholds
% ----------------------------------------------------
dur = sacc_t_ends - sacc_t_starts;

remove_idx = find( dur < options.min_sacc_duration | ...
                   dur > options.max_sacc_duration | ...
                   amplitudes < options.amplitude_thres );

if ~isempty(remove_idx)
    sacc_start_idx(remove_idx) = [];
    sacc_end_idx(remove_idx)   = [];
    sacc_t_starts(remove_idx)  = [];
    sacc_t_ends(remove_idx)    = [];
    amplitudes(remove_idx)     = [];
    dur(remove_idx)            = [];
end

durations = dur;


%% ----------------------------------------------------
%  Compute mean and peak velocities
% ----------------------------------------------------
peak_velocity = arrayfun(@(i) ...
    max(velocity(sacc_start_idx(i):sacc_end_idx(i))), ...
    1:length(sacc_end_idx));

mean_velocity = arrayfun(@(i) ...
    mean(velocity(sacc_start_idx(i):sacc_end_idx(i)), 'omitnan'), ...
    1:length(sacc_end_idx));


end  % function end
