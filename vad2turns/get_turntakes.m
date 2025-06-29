function [t,vspstate,tovl] = get_turntakes( varargin )
% get_turntakes - filter and classify turn takes
%
% [t,vspstate,tovl] = get_turntakes( t_start1, t_end1, t_start2, t_end2, ... )
%
% t_start : start times of speaking intervals
% t_end : end times of speaking intervals
%
% t :
% vspstate : sparse matrix of speaker activity with time stamps as first column
% tovl :
  vspstate = [];
  Nspeaker = numel(varargin)/2;
  % collect times of all interlocutors:
  for k_interlocutor=1:Nspeaker
    tstart = varargin{2*k_interlocutor-1};
    tend = varargin{2*k_interlocutor};
    % remove gaps as suggested by Eline Petersen 2023:
    [tstart, tend] = remove_gaps( tstart, tend, 1 );
    for k=1:numel(tstart)
      state = -ones(2,Nspeaker);
      state(1,k_interlocutor) = 1;
      state(2,k_interlocutor) = 0;
      state = [[tstart(k);tend(k)],state];
      vspstate(end+(1:2),1:(Nspeaker+1)) = state;
    end
  end
  % order sections:
  [tmp,idx]= sort(vspstate(:,1));
  vspstate = vspstate(idx,:);
  vspstate = [zeros(1,1+Nspeaker);vspstate;zeros(1,1+Nspeaker)];
  for k=2:size(vspstate,1)
    for kl=1:Nspeaker
      if vspstate(k,1+kl) < 0
        vspstate(k,1+kl) = vspstate(k-1,1+kl);
      end
    end
  end
  % find overlaps:
  idx_ovl = find(sum(vspstate(:,2:end),2)>1);
  tovl = [vspstate(idx_ovl,1),vspstate(idx_ovl+1,1)];
  % remove full overlaps (as in Petersen):
  idx_fullovl = idx_ovl(any(vspstate(idx_ovl,2:end) & (1-vspstate(idx_ovl-1,2:end)) & (1-vspstate(idx_ovl+1,2:end)),2));
  vspstateovl = vspstate(unique([idx_fullovl;idx_fullovl+1]),:);
  vspstate(unique([idx_fullovl;idx_fullovl+1]),:) = [];
  % collect overlap take time stamps:
  idx = find(any(vspstateovl(1:end-1,2:end) & (1-vspstateovl(2:end,2:end)),2));
  spkmatovl = double(vspstateovl(idx,2:end) & (1-vspstateovl(idx+1,2:end)));
  for k=2:size(spkmatovl,2)
    spkmatovl(:,k) = spkmatovl(:,k)*k;
  end
  tovl = [vspstateovl(idx,1),vspstateovl(idx+1,1),sum(spkmatovl,2)];
  % collect turn take time stamps:
  idx = find(any(vspstate(2:end,2:end) & (1-vspstate(1:end-1,2:end)),2))+1;
  spkmat = double(vspstate(idx,2:end) & (1-vspstate(idx-1,2:end)));
  for k=2:size(spkmat,2)
    spkmat(:,k) = spkmat(:,k)*k;
  end
  t = [vspstate(idx,1),vspstate(idx+1,1),sum(spkmat,2)];
  for k=1:size(t,1)
    speaker = t(k,3);
    idx0 = find(vspstate(idx(k):end,speaker+1)==0,1);
    t(k,2) = vspstate(idx(k)+idx0-1,1);
  end
  %vspstate
end
