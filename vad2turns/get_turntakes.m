function [t, mSpeakerState, tovl] = get_turntakes(varargin)
% GET_TURNTAKES Classify and filter turn-taking events from speaker intervals.
%
%   [t, mSpeakerState, tovl] = get_turntakes(t_start1, t_end1, t_start2, t_end2, ...)
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
%       mSpeakerState - Sparse matrix indicating speaker activity over time
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
%       [t, mSpeakerState, tovl] = get_turntakes(t_start1, t_end1, t_start2, t_end2);
%
%   See also: vad2turns, remove_gaps

% Check if the number of input arguments is even
    if mod(numel(varargin),2) ~= 0
        error('get_turntakes:OddNumberArguments', ...
              'An even number of arguments is required. Each speaker must have both start and end times.');
    end

    % Initialize output variables
    mSpeakerState = [];
    % Calculate the number of speakers based on input arguments
    Nspeaker = numel(varargin) / 2;

    % Collect and process speaking intervals for each speaker
    for k_interlocutor = 1:Nspeaker
        % Extract start and end times for the current speaker
        tstart = varargin{2 * k_interlocutor - 1};
        tend = varargin{2 * k_interlocutor};

        % Remove gaps of up to one second as suggested by Eline Petersen 2023
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
            % Append to mSpeakerState
            mSpeakerState(end + (1:2), 1:(Nspeaker + 1)) = state;
        end
    end

    % Sort the state sections by time
    [tmp, idx] = sort(mSpeakerState(:, 1));
    mSpeakerState = mSpeakerState(idx, :);

    % Add boundaries at the beginning and end
    mSpeakerState = [zeros(1, 1 + Nspeaker); mSpeakerState; zeros(1, 1 + Nspeaker)];

    % Fill in missing states between sections
    for k = 2:size(mSpeakerState, 1)
        for kl = 1:Nspeaker
            if mSpeakerState(k, 1 + kl) < 0
                mSpeakerState(k, 1 + kl) = mSpeakerState(k - 1, 1 + kl);
            end
        end
    end

    % Identify overlapping speech periods, i.e., any interval in which
    % more than one speaker is active
    idx_ovl = find(sum(mSpeakerState(:, 2:end), 2) > 1);
    tovl = [mSpeakerState(idx_ovl, 1), mSpeakerState(idx_ovl + 1, 1)];
    
    % Remove full overlaps (as in Petersen's method)
    idx_fullovl = [];
    for k=idx_ovl'
        is_ovl = true;
        % someone is speaking already:
        is_ovl = is_ovl & any(mSpeakerState(k,2:end) & mSpeakerState(k-1,2:end));
        % the one(s) who are speaking before are also speaking after:
        is_ovl = is_ovl & any(mSpeakerState(k-1,2:end) & mSpeakerState(k+1,2:end));
        if( is_ovl )
            idx_fullovl(end+1) = k;
        end
    end
    %idx_fullovl = idx_ovl(any(mSpeakerState(idx_ovl, 2:end) & (1 - mSpeakerState(idx_ovl - 1, 2:end)) & (1 - mSpeakerState(idx_ovl + 1, 2:end)), 2));
    mSpeakerStateovl = mSpeakerState(unique([idx_fullovl; idx_fullovl + 1]), :);
    mSpeakerState(unique([idx_fullovl; idx_fullovl + 1]), :) = [];

    % Collect overlap take time stamps
    idx = find(any(mSpeakerStateovl(1:end - 1, 2:end) & (1 - mSpeakerStateovl(2:end, 2:end)), 2));
    spkmatovl = double(mSpeakerStateovl(idx, 2:end) & (1 - mSpeakerStateovl(idx + 1, 2:end)));
    for k = 2:size(spkmatovl, 2)
        spkmatovl(:, k) = spkmatovl(:, k) * k;
    end
    tovl = [mSpeakerStateovl(idx, 1), mSpeakerStateovl(idx + 1, 1), sum(spkmatovl, 2)];

    % Collect turn take time stamps
    idx = find(any(mSpeakerState(2:end, 2:end) & (1 - mSpeakerState(1:end - 1, 2:end)), 2)) + 1;
    spkmat = double(mSpeakerState(idx, 2:end) & (1 - mSpeakerState(idx - 1, 2:end)));
    for k = 2:size(spkmat, 2)
        spkmat(:, k) = spkmat(:, k) * k;
    end
    t = [mSpeakerState(idx, 1), mSpeakerState(idx + 1, 1), sum(spkmat, 2)];

    % Adjust end times based on subsequent speaker activity
    for k = 1:size(t, 1)
        speaker = t(k, 3);
        idx0 = find(mSpeakerState(idx(k):end, speaker + 1) == 0, 1);
        t(k, 2) = mSpeakerState(idx(k) + idx0 - 1, 1);
    end
end
