function [tstart,tend] = remove_short_segments( tstart, tend, dmin )
  if nargin < 3
    dmin = 0.09;
    % value from Heldner2010
  end
  % length filtering:
  dur = tend-tstart;
  idx = find(dur >= dmin);
  tstart = tstart(idx);
  tend = tend(idx);
end