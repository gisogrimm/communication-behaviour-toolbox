function [gaze_in_roi, sortedMeans] = gaze_roi_gmm(gaze_az, nROI)
% Identify gaze ROIs using GMM 
%
% INPUTS:
%   gaze_az : nSamples × 1 vector of horizontal gaze positions (azimuth)
%   nROI    : starting number of Gaussians (e.g., 3)
%
% OUTPUT:
%   gaze_in_roi : nSamples × nroi binary matrix
%                 1 if the gaze sample belongs to that ROI, 0 otherwise
%
% NOTES:
%   - The function fits a GMM starting with nROI components and reduces if
%     it doesn't converge.
%   - for nroi=3 Columns are assigned according to gaze azimuth:
%       * Col1 = leftmost
%       * Col2 = rightmost
%       * Col3 = middle
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

maxTry = nROI;

%% Fit GMM, decreasing components if it fails
while maxTry > 0
    try
        gm = fitgmdist(gaze_az, maxTry, 'Options', statset('MaxIter',1000,'Display','off'));
        if gm.Converged
            break;
        else
            maxTry = maxTry - 1;
        end
    catch
        maxTry = maxTry - 1;
    end
end

%% If no GMM converged, return NaNs
if maxTry == 0
    gaze_in_roi = NaN(length(gaze_az),3);
    return;
end

%% Predict cluster assignments
idx = cluster(gm, gaze_az);
nFitted = max(idx);  % actual number of clusters

% Cluster means
cluster_means = gm.mu;
    [sortedMeans, sortIdx] = sort(cluster_means);

%% Convert to binary matrix: nSamples × nFitted
if nROI==3

%% Initialize 3-column output
gaze_in_roi = zeros(length(gaze_az), 3);

%% Assign clusters to columns
if nFitted == 3
    % Sort means: left < middle < right
    col1 = sortIdx(1); % leftmost
    col3 = sortIdx(2); % middle
    col2 = sortIdx(3); % rightmost
    
    gaze_in_roi(:,1) = idx == col1;
    gaze_in_roi(:,2) = idx == col2;
    gaze_in_roi(:,3) = idx == col3;

elseif nFitted == 2
    col1 = sortIdx(1); % left
    col2 = sortIdx(2); % right
    % Assign first two columns, third remains zeros
    gaze_in_roi(:,1) = idx == col1;
    gaze_in_roi(:,2) = idx == col2;
    gaze_in_roi(:,3) = 0;

elseif nFitted == 1
    sortedMeans = cluster_means(1);
    if sortedMeans < 0
        gaze_in_roi(:,2) = idx == 1; % negative → col2
    else
        gaze_in_roi(:,1) = idx == 1; % positive → col1
    end
    gaze_in_roi(:,3) = 0;
end
else
    gaze_in_roi = zeros(length(gaze_az), nFitted);

for r = 1:nFitted
    gaze_in_roi(idx == r, r) = 1;
end
end

end
