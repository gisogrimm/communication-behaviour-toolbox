%NB this function used in mergeAndRemoveBurst.m
function [on, off] = detectChanges(actArr)

if size(actArr,2)< size(actArr,1)
    actArr=actArr';
end

% Determine sign changes
signChange = diff([actArr(1) actArr]); % +1 on, -1: off

% Time indices of activity on and off
on  = find(signChange == 1);
off = find(signChange == -1);

% Define first and last sign change
[~,~,firstSgnCgn] = find(signChange,1,'first');
[~,~,lastSgnCgn]  = find(signChange,1,'last');

% Add "on" signal to beginning of array
if firstSgnCgn == -1
    on = [1 on];
end

% Add "off" signal to end of array
if lastSgnCgn == 1
    off = [off length(actArr)];
end

end

