function [talking, listening] = downsampleSpeech(gazeTime, turns, talkerIdx)
% DOWNSAMPLESPEECH  Create talking/listening vectors for gaze timestamps
%
% INPUTS:
%   gazeTime  : nSamples × 1 vector of time stamps
%   turns     : either
%               - struct array with turns(i).turns = [nTurns x 2] start/end times
%               - matrix [nTurns x 3] = [start end talkerIdx]
%   talkerIdx : integer index of the talker of interest
%
% OUTPUTS:
%   talking   : nSamples × 1 binary vector → 1 if talkerIdx is talking
%   listening : nSamples × nOtherTalkers binary matrix → 1 if talkerIdx is listening to the other talker
%

%% Normalize gazeTime
gazeTime = gazeTime(:);
nSamples = length(gazeTime);

%% Convert turns to struct if needed
if ~isstruct(turns)
    nTalkers = length(unique(turns(:,3)));
    tmpTurns = turns;
    turnsStruct = struct();
    for t = 1:nTalkers
        turnsStruct(t).turns = tmpTurns(tmpTurns(:,3) == t, 1:2);
    end
    turns = turnsStruct;
else
    nTalkers = length(turns);
end

%% Initialize output
talking = zeros(nSamples,1);

% All other talkers for listening
otherTalkers = setdiff(1:nTalkers, talkerIdx);
nOther = length(otherTalkers);
listening = zeros(nSamples, nOther);

%% Fill talking vector
myTurns = turns(talkerIdx).turns;  % [start end] for this talker
for k = 1:size(myTurns,1)
    idx = gazeTime >= myTurns(k,1) & gazeTime <= myTurns(k,2);
    talking(idx) = 1;
end

%% Fill listening matrix
for j = 1:nOther
    tIdx = otherTalkers(j);
    otherTurns = turns(tIdx).turns;
    for k = 1:size(otherTurns,1)
        idx = gazeTime >= otherTurns(k,1) & gazeTime <= otherTurns(k,2);
        listening(idx,j) = 1;
    end
end

end
