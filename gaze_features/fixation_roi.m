function fix_in_roi = fixation_roi(fix, gaze_in_roi)
% Determines whether each fixation lies within a ROI
%
% INPUTS:
%   fix          : n×1 vector of fixation labels per gaze sample
%                  0 or NaN = not part of a fixation
%                  1,2,3,... = fixation index
%
%   gaze_in_roi  : n×1 logical or numeric vector
%                  1 = gaze sample inside ROI
%                  0 = outside ROI
%
% OUTPUT:
%   fix_in_roi   : f×1 logical vector
%                  f = number of fixations
%                  1 = fixation is in ROI
%                  0 = fixation is outside ROI
%
% A fixation is considered "in ROI" if the majority of its samples
% fall within the ROI.
%
arguments
    fix double
    gaze_in_roi double
end

%% Input normalization
fix = fix(:);
gaze_in_roi = gaze_in_roi(:);

if length(fix) ~= length(gaze_in_roi)
    error('fix and gaze_in_roi must be the same length');
end

%% Get unique fixation indices
fix_ids = unique(fix);
fix_ids(isnan(fix_ids) | fix_ids == 0) = [];  % remove invalid labels

n_fix = length(fix_ids);
fix_in_roi = false(n_fix,1);

%% Determine ROI membership per fixation
for i = 1:n_fix
    idx = fix == fix_ids(i);            % samples belonging to fixation i
    
    % if any then yes
    fix_in_roi(i) = any(gaze_in_roi(idx));
end

end
