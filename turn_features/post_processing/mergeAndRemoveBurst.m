function [actArr,idx] = mergeAndRemoveBurst(actArr,timeVec,burstThresh,silThresh)
%merge speech segments with small gaps and remove bursts
arguments
    actArr logical %speech activity vector
    timeVec double %time vector
    burstThresh double = 90e-3;                                                % Threshold for small bursts that are unlikely to be speech
    silThresh   double = 180e-3;                                               % Threshold for defining silent gaps
end


% Set first sample to zero to detect changes in case it's active from the
% beginning:
actArr(1) = 0;
actArr(end) = 0;


% BRIDGE GAPS:
%
[on, off] = detectChanges(actArr);
%
% Define gap duration
gapDuration = abs(timeVec(off(1:end-1))-timeVec(on(2:end)));
%
% Remove pauses < 180 ms
for j = 1:length(gapDuration)
    if gapDuration(j) <= silThresh
        actArr(off(j):on(j+1)) = 1;
    end
end


% REMOVE IRRELEVANT BURSTS OF SPEECH:
%
% Detect changes after bursts within a stream are bridged
%
[on, off] = detectChanges(actArr);
%
% Define burst duration
burst = timeVec(off)-timeVec(on);
%
% Set bursts < 90 ms to zero
for i = 1:length(burst)
    if burst(i) <= burstThresh
        actArr(on(i):off(i)) = 0;

    end
end

idx=[find(actArr(2:end)-actArr(1:end-1)==1);find(actArr(2:end)-actArr(1:end-1)==-1)]';
end