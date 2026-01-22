function roi_out = postprocessROI(roi, t, minDur, minGap)
% Temporal processing of ROI gaze vector
%
% roi_out = processROI(roi, t, minDur, minGap)
%
% roi     : vector of ROI indices (0 = no ROI)
% t       : time vector (same length as roi)
% minDur  : minimum duration of an ROI segment
% minGap  : maximum gap duration to merge ROI segments
%
% roi_out : processed ROI vector
arguments
    roi double
    t double 
    minDur double = 0.1 %seconds
    minGap double = 0.1 %seconds
end
roi = roi(:);
t   = t(:);

roi_out = roi;
N = length(roi);

%% -------------------------------------------------
% 1) Merge short gaps between same ROI
%% -------------------------------------------------
i = 1;
while i < N
    if roi_out(i) > 0
        r = roi_out(i);

        % find end of this ROI segment
        i_end = i;
        while i_end < N && roi_out(i_end+1) == r
            i_end = i_end + 1;
        end

        % look ahead for same ROI after a gap
        j = i_end + 1;
        while j < N && roi_out(j) == 0
            j = j + 1;
        end

        if j <= N && roi_out(j) == r
            gapDur = t(j) - t(i_end);
            if gapDur <= minGap
                % fill the gap
                roi_out(i_end+1:j-1) = r;
            end
        end

        i = i_end + 1;
    else
        i = i + 1;
    end
end

%% -------------------------------------------------
% 2) Remove short ROI segments
%% -------------------------------------------------
i = 1;
while i <= N
    if roi_out(i) > 0
        r = roi_out(i);

        % find end of this ROI segment
        i_end = i;
        while i_end < N && roi_out(i_end+1) == r
            i_end = i_end + 1;
        end

        dur = t(i_end) - t(i);
        if dur < minDur
            roi_out(i:i_end) = 0;
        end

        i = i_end + 1;
    else
        i = i + 1;
    end
end
end
