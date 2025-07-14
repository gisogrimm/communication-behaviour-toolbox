function [t, vspstate, tovl] = get_turntakes(varargin)
% GET_TURNTAKES Classify and filter turn-taking events from speaker intervals.
%
%   [t, vspstate, tovl] = get_turntakes(t_start1, t_end1, t_start2, t_end2, ...)
%   processes the start and end times of speaking intervals for multiple
%   speakers to classify turn-taking events and overlaps.
%
%   Inputs:
%       t_start1, t_end1, ... - Pairs of start and end times for each speaker
%       (optional arguments, must be even number of arguments)
%
%   Outputs:
%       t - Matrix of turn-taking events with columns:
%           [start_time, end_time, speaker_id]
%       vspstate - Sparse matrix indicating speaker activity over time
%       tovl - Matrix of overlap events with columns:
%           [start_time, end_time, overlapping_speakers]
%
%   The function processes the input intervals to:
%   - Remove gaps in speaking intervals
%   - Classify turn-taking events
%   - Identify overlapping speech periods
%
%   Example:
%       % Sample usage with two speakers
%       t_start1 = [0, 2]; t_end1 = [1, 3];
%       t_start2 = [0.5, 1.5]; t_end2 = [2, 3.5];
%       [t, vspstate, tovl] = get_turntakes(t_start1, t_end1, t_start2, t_end2);
%
%   See also: vad2turns, remove_gaps

% Check if the number of input arguments is even
    if mod(numel(varargin),2) ~= 0
        error('get_turntakes:OddNumberArguments', ...
              'An even number of arguments is required. Each speaker must have both start and end times.');
    end

    % Initialize output variables
    vspstate = [];
    % Calculate the number of speakers based on input arguments
    Nspeaker = numel(varargin) / 2;

    % Collect and process speaking intervals for each speaker
    for k_interlocutor = 1:Nspeaker
        % Extract start and end times for the current speaker
        tstart = varargin{2 * k_interlocutor - 1};
        tend = varargin{2 * k_interlocutor};

        % Remove gaps as suggested by Eline Petersen 2023
        [tstart, tend] = remove_gaps(tstart, tend, 1);

        % Process each speaking interval
        for k = 1:numel(tstart)
            % Initialize state matrix with -1s (inactive)
            state = -ones(2, Nspeaker);
            % Set the current speaker as active (1) and inactive (0) at end
            state(1, k_interlocutor) = 1;
            state(2, k_interlocutor) = 0;
            % Combine time and state information
            state = [[tstart(k); tend(k)], state];
            % Append to vspstate
            vspstate(end + (1:2), 1:(Nspeaker + 1)) = state;
        end
    end

    % Sort the state sections by time
    [tmp, idx] = sort(vspstate(:, 1));
    vspstate = vspstate(idx, :);

    % Add boundaries at the beginning and end
    vspstate = [zeros(1, 1 + Nspeaker); vspstate; zeros(1, 1 + Nspeaker)];

    % Fill in missing states between sections
    for k = 2:size(vspstate, 1)
        for kl = 1:Nspeaker
            if vspstate(k, 1 + kl) < 0
                vspstate(k, 1 + kl) = vspstate(k - 1, 1 + kl);
            end
        end
    end

    % Identify overlapping speech periods
    idx_ovl = find(sum(vspstate(:, 2:end), 2) > 1);
    tovl = [vspstate(idx_ovl, 1), vspstate(idx_ovl + 1, 1)];

    % Remove full overlaps (as in Petersen's method)
    idx_fullovl = idx_ovl(any(vspstate(idx_ovl, 2:end) & (1 - vspstate(idx_ovl - 1, 2:end)) & (1 - vspstate(idx_ovl + 1, 2:end)), 2));
    vspstateovl = vspstate(unique([idx_fullovl; idx_fullovl + 1]), :);
    vspstate(unique([idx_fullovl; idx_fullovl + 1]), :) = [];

    % Collect overlap take time stamps
    idx = find(any(vspstateovl(1:end - 1, 2:end) & (1 - vspstateovl(2:end, 2:end)), 2));
    spkmatovl = double(vspstateovl(idx, 2:end) & (1 - vspstateovl(idx + 1, 2:end)));
    for k = 2:size(spkmatovl, 2)
        spkmatovl(:, k) = spkmatovl(:, k) * k;
    end
    tovl = [vspstateovl(idx, 1), vspstateovl(idx + 1, 1), sum(spkmatovl, 2)];

    % Collect turn take time stamps
    idx = find(any(vspstate(2:end, 2:end) & (1 - vspstate(1:end - 1, 2:end)), 2)) + 1;
    spkmat = double(vspstate(idx, 2:end) & (1 - vspstate(idx - 1, 2:end)));
    for k = 2:size(spkmat, 2)
        spkmat(:, k) = spkmat(:, k) * k;
    end
    t = [vspstate(idx, 1), vspstate(idx + 1, 1), sum(spkmat, 2)];

    % Adjust end times based on subsequent speaker activity
    for k = 1:size(t, 1)
        speaker = t(k, 3);
        idx0 = find(vspstate(idx(k):end, speaker + 1) == 0, 1);
        t(k, 2) = vspstate(idx(k) + idx0 - 1, 1);
    end
end
