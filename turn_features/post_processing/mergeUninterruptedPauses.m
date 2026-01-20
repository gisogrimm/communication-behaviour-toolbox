function [act, idx] = mergeUninterruptedPauses(act,idx,thres)
%function to merge uninterrupted pauses shorter than thres
arguments
    act logical 
    idx struct
    thres double
end
   

if nargin<3
    thres=100;
end

if size(act,1)>size(act,2)
    act=act';
end

talkers=1:size(act,1);
for subj = 1:size(act,1)
    [on, off] = detectChanges(act(subj,:));

    %check if someone else is talking
    for i=1:length(off)-1
        interl=talkers(talkers~=subj);%other talker indices
        %check if someone else is talking in pause
        if ~any(act(interl,off(i):on(i+1)))
            %check if pause is less than break
            if on(i+1)-off(i)<thres
            act(subj,off(i):on(i+1))=1;
            end
        end
    end

    %return index struct as well
    if nargout > 1
        idx(subj).turns=[find(act(subj,2:end)-act(subj,1:end-1)==1);find(act(subj,2:end)-act(subj,1:end-1)==-1)]';
    end
end


end