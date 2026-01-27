function turns = computeTurnFeatures(turns, VADarr)
%% COMPUTETURNFEATURES
% This function computes turn-level conversational features:
%
%   - Previous talker and time-from-previous-offset (FTOstart)
%   - Next talker and time-to-next-onset (FTOend)
%   - Turn duration
%
% INPUTS:
%   turns  : Either
%            (1) struct array, one element per talker, with field:
%                   turns(i).turns = [nTurns × 2] [start end] (sample indices)
%            OR
%            (2) numeric matrix [nTurns × 3]:
%                   col 1 = turn start
%                   col 2 = turn end
%                   col 3 = talker index
%
%   VADarr : [nTalkers × nSamples] binary voice activity matrix
%            VADarr(talker, sample) == 1 → talker is active
%
% OUTPUT:
%   turns  : Same format as input, enriched with:
%            - prevTalker
%            - nextTalker
%            - FTOstart (time from previous offset)
%            - FTOend   (time to next onset)
%            - turnDur
%

%% ------------------------------------------------------------
% Check input format and convert matrix → struct if needed
% ------------------------------------------------------------
transform = ~isstruct(turns); % remember whether we must convert back

if nargin<2
    error('not enought input arguments')
end

if transform
    % Number of unique talkers
    nTalker = length(unique(turns(:,3)));

    % Backup original matrix
    tmpTurns = turns;

    % Convert matrix format to struct array
    turns = struct();
    for i = nTalker
        turns(i).turns = tmpTurns(tmpTurns(:,3) == i, 1:2);
    end
end

if size(VADarr,1)>size(VADarr,2)
    VADarr=VADarr';
end

%% ------------------------------------------------------------
% Compute NEXT TALKER and FTOend for each turn
% ------------------------------------------------------------
for i = 1:length(turns)

    % Preallocate output fields
    turns(i).nextTalker = zeros(size(turns(i).turns,1),1);

    for t = 1:size(turns(i).turns,1)

        % tmp(j) = distance (in samples) until talker j speaks next
        tmp = zeros(length(turns),1);

        for j = 1:length(turns)
            try
                tmp(j) = find(VADarr(j, turns(i).turns(t,2)+1:end) == 1, 1);
            catch
                tmp(j) = NaN;
            end
        end

        % No following speech detected
        if all(isnan(tmp))
            turns(i).nextTalker(t) = NaN;
            turns(i).FTOend(t) = NaN;

        else
            % Find earliest next speaker
            ch = find(tmp == min(tmp));

            % Resolve ties
            if length(ch) > 1
                for c = ch'
                    a = find( ...
                        VADarr(c, 1:turns(i).turns(t,2)+1+min(tmp)) == 0, ...
                        1, 'last') + 1;

                    if ~isempty(turns(c).turns)
                        if ismember(a, turns(c).turns(:,1))
                            turns(i).nextTalker(t) = c;
                            turns(i).FTOend(t) = a - turns(i).turns(t,2);
                        end
                    end
                end

                % If unresolved
                if turns(i).nextTalker(t) == 0
                    turns(i).nextTalker(t) = NaN;
                    turns(i).FTOend(t) = NaN;
                end

            else
                % Unique next speaker
                [minIdx, turns(i).nextTalker(t)] = min(tmp);
                turns(i).FTOend(t) = ...
                    find(VADarr(turns(i).nextTalker(t), ...
                    1:turns(i).turns(t,2)+1+minIdx) == 0, 1, 'last') ...
                    + 1 - turns(i).turns(t,2);
            end
        end
    end

    %% --------------------------------------------------------
    % Compute PREVIOUS TALKER and FTOstart for each turn
    % --------------------------------------------------------
    turns(i).prevTalker = zeros(size(turns(i).turns,1),1);

    for t = 1:size(turns(i).turns,1)

        % tmp(j) = last time talker j spoke before current turn
        tmp = zeros(length(turns),1);

        for j = 1:length(turns)
            try
                tmp(j) = find(VADarr(j, 1:turns(i).turns(t,1)-1) == 1, 1, 'last');
            catch
                tmp(j) = NaN;
            end
        end

        % No previous speaker detected
        if all(isnan(tmp))
            turns(i).prevTalker(t) = NaN;
            turns(i).FTOstart(t) = NaN;

        else
            % Find most recent previous speaker
            ch = find(tmp == max(tmp));

            % Resolve ties
            if length(ch) > 1
                for c = ch'
                    a = find( ...
                        VADarr(c, max(tmp):turns(i).turns(t,2)) == 0, ...
                        1) + max(tmp) - 1;

                    if ~isempty(turns(c).turns)
                        if ismember(a, turns(c).turns(:,2))
                            turns(i).prevTalker(t) = c;
                            turns(i).FTOstart(t) = turns(i).turns(t,1) - a;
                        end
                    end
                end

                % If unresolved
                if turns(i).prevTalker(t) == 0
                    turns(i).prevTalker(t) = NaN;
                    turns(i).FTOstart(t) = NaN;
                end

            else
                % Unique previous speaker
                [maxIdx, turns(i).prevTalker(t)] = max(tmp);
                turns(i).FTOstart(t) = ...
                    turns(i).turns(t,1) - ...
                    find(VADarr(turns(i).prevTalker(t), ...
                    maxIdx:turns(i).turns(t,2)) == 0, 1) - maxIdx + 1;
            end
        end
    end

    %% --------------------------------------------------------
    % Turn duration
    % --------------------------------------------------------
    turns(i).turnDur = turns(i).turns(:,2) - turns(i).turns(:,1);
end

%% ------------------------------------------------------------
% Convert back to matrix format if needed
% ------------------------------------------------------------
if transform
    tmpTurns(:,4) = NaN; % prevTalker
    tmpTurns(:,5) = NaN; % FTOstart
    tmpTurns(:,6) = NaN; % nextTalker
    tmpTurns(:,7) = NaN; % FTOend

    for i = nTalker
        tmpTurns(tmpTurns(:,3) == i, 4) = turns(i).prevTalker;
        tmpTurns(tmpTurns(:,3) == i, 5) = turns(i).FTOstart;
        tmpTurns(tmpTurns(:,3) == i, 6) = turns(i).nextTalker;
        tmpTurns(tmpTurns(:,3) == i, 7) = turns(i).FTOend;
    end

    turns = tmpTurns;
end
end
