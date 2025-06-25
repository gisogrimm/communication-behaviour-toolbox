function [tstart, tend] = remove_gaps( tstart, tend, dmin )
  if nargin < 3
    dmin = 0.18;
    % value from Heldner2010
  end
  % gap filtering:
  gap_dur = tstart(2:end)-tend(1:end-1);
  idx = find(gap_dur < dmin);
  tstart(idx+1) = [];
  tend(idx) = [];
end