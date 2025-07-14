function [mVADout, t, vspstate, tovl] = vad2turns(vT, mVAD, maxgapdur, minsegmentdur)
% VAD2TURNS Processes voice activity detection data to identify speaker turns.
%
%   [mVADout, t, vspstate, tovl] = vad2turns(vT, mVAD, maxgapdur, minsegmentdur)
%   processes the voice activity detection (VAD) data in mVAD to identify
%   speaker turns based on the timing information in vT. The function
%   accounts for maximum gap duration (maxgapdur) and minimum segment
%   duration (minsegmentdur) to refine the turn detection.
%
%   Inputs:
%       vT - Time vector corresponding to the VAD data (vector)
%       mVAD - Voice activity detection matrix (matrix)
%           - Rows represent time points
%           - Columns represent different speakers
%       maxgapdur - Maximum allowable gap duration between speech segments
%                   (default: 0.18 seconds, from Heldner et al. 2010)
%       minsegmentdur - Minimum duration of a speech segment to be considered
%                       (default: 0.09 seconds, from Heldner et al. 2010)
%
%   Outputs:
%       mVADout - Refined VAD matrix with updated speech segments
%       t - Time points of speaker turns
%       vspstate - Speaker state information
%       tovl - Overlap information between speaker turns
%
%   The function processes each channel to detect active speech periods,
%   removes gaps shorter than maxgapdur, filters out short segments
%   shorter than minsegmentdur, and updates the VAD matrix accordingly.
%   It then determines turn-taking information based on the refined segments.
%
%   Example:
%       % Sample usage with default parameters
%       [mVADout, t, vspstate, tovl] = vad2turns(vT, mVAD);
%
%   See also remove_gaps, remove_short_segments, get_turntakes

% Set default values for maxgapdur and minsegmentdur if not provided
    if nargin < 3
        maxgapdur = 0.18; % Default from Heldner et al. 2010
    end
    if nargin < 4
        minsegmentdur = 0.09; % Default from Heldner et al. 2010
    end

    % Initialize the output VAD matrix with zeros
    mVADout = zeros(size(mVAD));

    % Prepare cell array to collect arguments for get_turntakes
    cArg = {};

    % Process each channel to detect and refine speech segments
    for ch = 1:size(mVAD, 2)
        % Extract activity for the current channel
        act = mVAD(:, ch);
        
        % Ensure the first and last elements are inactive to handle edge cases
        act(1) = 0;
        act(end + 1) = 0;
        
        % Calculate differences to find transitions
        dact = diff(act);
        
        % Find start and end times of active segments
        tstart = vT(find(dact == 1)); % Start times where activity begins
        tend = vT(find(dact == -1));  % End times where activity ceases
        
        % Remove gaps shorter than maxgapdur
        [tstart, tend] = remove_gaps(tstart, tend, maxgapdur);
        
        % Filter out short segments based on minsegmentdur
        [tstart, tend] = remove_short_segments(tstart, tend, minsegmentdur);
        
        % Update mVADout for the current channel
        for k = 1:numel(tstart)
            % Find indices corresponding to start and end times
            idx_start = find(vT == tstart(k), 1);
            idx_end = find(vT == tend(k), 1);
            
            % Set the VAD to active (1) for the duration of the segment
            mVADout(idx_start:idx_end, ch) = 1;
        end
        
        % Collect start and end times for turn-taking analysis
        cArg{end + 1} = tstart;
        cArg{end + 1} = tend;
    end

    % Determine turn-taking information
    [t, vspstate, tovl] = get_turntakes(cArg{:});
end
