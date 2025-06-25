vT = readmatrix('group1/quiet/T.csv');
mVAD = readmatrix('group1/quiet/VAD.csv');
[mVADsmooth,t_turns,vspstate,tovl] = vad2turns( vT, mVAD);


figure
mCol = colororder();
for k=1:size(t_turns,1)
    spk = t_turns(k,3);
    patch([t_turns(k,1),t_turns(k,2),t_turns(k,2),t_turns(k,1)],...
          [1,1,0,0]-0.5+spk,mCol(spk,:));
end
set(gca,'YTick',1:size(mVAD,2));
