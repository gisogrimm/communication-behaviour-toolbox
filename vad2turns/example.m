% read the time vector:
vT = readmatrix('group1/quiet/T.csv');
% ead voice activity matrix:
mVAD = readmatrix('group1/quiet/VAD.csv');
% convert voice activity into the sparse matrix t_turns:
[mVADsmooth,t_turns,vspstate,tovl] = vad2turns( vT, mVAD);

% plot turns:
figure
mCol = colororder();
for k=1:size(t_turns,1)
    spk = t_turns(k,3);
    patch([t_turns(k,1),t_turns(k,2),t_turns(k,2),t_turns(k,1)],...
          [1,1,0,0]-0.5+spk,mCol(spk,:));
end
set(gca,'YTick',1:size(mVAD,2));
