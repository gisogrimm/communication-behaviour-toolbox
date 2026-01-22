function [gaze_in_roi, sortedMeans] = gaze_roi_kmeans(gaze_az, nROI)
% Identify gaze ROIs using GMM
%
% INPUTS:
%   gaze_az : nSamples × 1 vector of horizontal gaze positions (azimuth)
%   nROI    : starting number of Gaussians (e.g., 3)
%
% OUTPUT:
%   gaze_in_roi : nSamples × 1  matrix
%                 one number for each ROI, 0 otherwise
%
% NOTES:
%   - The function fits a GMM starting with nROI components and reduces if
%     it doesn't converge.
%   - for nroi=3 Columns are assigned according to gaze azimuth:
%       * 1 = leftmost
%       * 2 = rightmost
%       * 3 = middle
%
arguments
    gaze_az double
    nROI double = 3
end

%% Input normalization
gaze_az = gaze_az(:);
if isempty(gaze_az)
    gaze_in_roi = NaN(1,3);
    return;
end

[idx,cluster_means]=kmeans(gaze_az,nROI);

%% Order
[sortedMeans, sortIdx] = sort(cluster_means);

if nROI==3

    %% Initialize 3-column output
    gaze_in_roi = zeros(length(gaze_az), 1);

    %% Assign clusters to columns
        % Sort means: left < middle < right
        gaze_in_roi( idx==sortIdx(1)) = 1;% leftmost
        gaze_in_roi( idx==sortIdx(2)) = 3;% middle
        gaze_in_roi(idx== sortIdx(3)) = 2; % rightmost
   
else
    gaze_in_roi = zeros(length(gaze_az), 1);

    for r = 1:nROI
        gaze_in_roi(idx == r) = r;
    end
end

end
