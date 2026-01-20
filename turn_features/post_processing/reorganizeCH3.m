function [CH_reorg,commClassCorrection] = reorganizeCH3(CH)
%% Reorganize turn information for triads
% based on output of this Conversational dynamic classifier
% Anna Josefine Munch Sørensen. (2025). Communicative State Classification: 
% MATLAB algorithm for classifying multi-talker audio into timings and 
% durations of conversational dynamics (Version 1.0.0) 
% [Computer software]. https://doi.org/10.5281/zenodo.17063065
%
% Correct the talker order of the CH and performs deletion of overlapping
% overlaps. The actions are as follows:
% 1) reorganize talkers such that overlaps occurs for the speaker in question
% 2) Correct for overlapping overlapB (overlapping two different talkers
% changed to overlapping the talker starting first)
% 3)Correct overlapping overlapB and overlapW for same talker --> remove overlapB
% 4) correct for overlapping overlapW --> remove the second entry (always
% identical in time and duration)
%
%5) IPU after (talkers own) gap and overlapB are calculated
%6) IPU before (other talkers) gap and overlapB are calculated
CH_org=CH;
commClassCorrection.multiOverlapB = 0;
commClassCorrection.multiOverlapW = 0;
commClassCorrection.multiOverlapBW = 0;
%% Reprganize classes to belong to the correct talker
classes2reorganize = {'gap';'overlapB';'overlapW'};

for i_class = 1:size(classes2reorganize,1)

    %init
    for i_talker =1:length(CH)
        oldClass{i_talker} = CH{i_talker}.(classes2reorganize{i_class});

        newClass{i_talker}.startIdx = [];
        newClass{i_talker}.endIdx = [];
        newClass{i_talker}.duration = [];
        newClass{i_talker}.channel = [];
        newClass{i_talker}.numData = [];
    end

    % reorganize
    for i_talker = 1:length(CH)
        othertalkers = 1:length(CH) ;
        othertalkers(i_talker)=[];

        for i_partner = 1:length(othertalkers)
            newClass{i_talker}.startIdx = [newClass{i_talker}.startIdx oldClass{othertalkers(i_partner)}.startIdx(oldClass{othertalkers(i_partner)}.channel == i_talker)];
            newClass{i_talker}.endIdx = [newClass{i_talker}.endIdx oldClass{othertalkers(i_partner)}.endIdx(oldClass{othertalkers(i_partner)}.channel == i_talker)];
            newClass{i_talker}.duration = [newClass{i_talker}.duration oldClass{othertalkers(i_partner)}.duration(oldClass{othertalkers(i_partner)}.channel == i_talker)];
            newClass{i_talker}.channel = [newClass{i_talker}.channel repmat(othertalkers(i_partner),1,sum(oldClass{othertalkers(i_partner)}.channel == i_talker))];
        end
    end

    % sort accordign to startIdx
    for i_talker = 1:length(CH)
        [~,sorting]=sort(newClass{i_talker}.startIdx);

        newClass{i_talker}.startIdx = newClass{i_talker}.startIdx(sorting);
        newClass{i_talker}.endIdx = newClass{i_talker}.endIdx(sorting);
        newClass{i_talker}.duration = newClass{i_talker}.duration(sorting);
        newClass{i_talker}.channel = newClass{i_talker}.channel(sorting);

        % if there for some reason is dublications of a gap remove it
        % if any(diff(newClass{i_talker}.startIdx)==0)
        %     tmp_rmindx = find(diff( newClass{i_talker}.startIdx)==0);
        % 
        %     newClass{i_talker}.startIdx(tmp_rmindx)=[];
        %     newClass{i_talker}.endIdx(tmp_rmindx)=[];
        %     newClass{i_talker}.duration(tmp_rmindx)=[];
        %     newClass{i_talker}.channel(tmp_rmindx)=[];
        % end
        newClass{i_talker}.numData = numel(newClass{i_talker}.startIdx);

        % put back into CH structure
        CH{i_talker}.(classes2reorganize{i_class})  = newClass{i_talker};

    end

end

%% 2) Correct for overlapping overlapBs
% correct multiple detections of overlapB for the same speaker --> extend
% overlap period to include both time intervals
% OBS: Some overlaps can other start times than the IPUs, like speaker 3
% below which will have two overlaps: 1) with talker 1 starting at IPU start
% and 2) with talker 2 with an IPU defined by the start of talker 2's IPU
%                       ---------------
%                               ---------------
%                            -----------
for i_talker =1:length(CH)
    % check if start of X and X+1 overlapB are shorter apart then the duration of X
    overlappingDetections = find(diff(CH{i_talker}.overlapB.startIdx) <= (CH{i_talker}.overlapB.duration(1:end-1)));
    overlaporg{i_talker}=CH{i_talker}.overlapB;
    % if overlapping overlapB's exist
    if ~isempty(overlappingDetections)

        commClassCorrection.multiOverlapB = length(overlappingDetections);

        % Delete one of the two overlaps and extend stop time if needed
        for i_overlap = 1:length(overlappingDetections)
            %             fprintf('B')
    
            
             overlaps=[CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)) CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)+1)];
 isTurn=find(ismember(overlaps,CH{i_talker}.turn.startIdx));
 if length(isTurn)>1
        %check end for the overlaps       
    overlapsE=[CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)) CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)+1)];
    overlapsC=[CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)) CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)+1)];
    isTurn=find([ismember(overlapsE(1),CH{overlapsC(1)}.turn.endIdx),ismember(overlapsE(2),CH{overlapsC(2)}.turn.endIdx)]);
 end

if ~isempty(isTurn) && length(isTurn)==1
            CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)) = overlaps(isTurn);
            CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)+(isTurn-1));
                  CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)+(isTurn-1));
            CH{i_talker}.overlapB.duration(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.duration(overlappingDetections(i_overlap)+(isTurn-1));
elseif length(isTurn)>1
 
    % start at min stop indx
            [tmp_minval,tmp_minindx] = ...
                min([CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)) CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)+1)]);
  CH{i_talker}.overlapB.startIdx(overlappingDetections(i_overlap)) = tmp_minval;
           CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.endIdx(overlappingDetections(i_overlap)+(tmp_minindx-1));
             CH{i_talker}.overlapB.duration(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.duration(overlappingDetections(i_overlap)+(tmp_minindx-1));
           CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)) = CH{i_talker}.overlapB.channel(overlappingDetections(i_overlap)+(tmp_minindx-1));
 
end
        end

        % delete overlapping overlap
        CH{i_talker}.overlapB.startIdx(overlappingDetections+1)=[];
        CH{i_talker}.overlapB.endIdx(overlappingDetections+1)=[];
        CH{i_talker}.overlapB.duration(overlappingDetections+1)=[];
        CH{i_talker}.overlapB.channel(overlappingDetections+1)=[];

        CH{i_talker}.overlapB.numData = numel(CH{i_talker}.overlapB.duration);

    end

    % check that overlap start is associated with an IPU start
    rmIdx=[];
    for i_overlap = 1:length(CH{i_talker}.overlapB.startIdx)
        if isempty(intersect(CH{i_talker}.overlapB.startIdx(i_overlap),CH{i_talker}.IPU.startIdx))
            rmIdx = [rmIdx i_overlap];
        end

        %and turn start
        if isempty(intersect(CH{i_talker}.overlapB.startIdx(i_overlap),CH{i_talker}.turn.startIdx))
            rmIdx = [rmIdx i_overlap];
        end
  %and turn start
  if isempty(intersect(CH{i_talker}.overlapB.endIdx(i_overlap),CH{CH{i_talker}.overlapB.channel(i_overlap)}.turn.endIdx))
            rmIdx = [rmIdx i_overlap];
        end
    end

    CH{i_talker}.overlapB.startIdx(rmIdx)=[];
    CH{i_talker}.overlapB.endIdx(rmIdx)=[];
    CH{i_talker}.overlapB.duration(rmIdx)=[];
    CH{i_talker}.overlapB.channel(rmIdx)=[];

    CH{i_talker}.overlapB.numData = numel(CH{i_talker}.overlapB.duration);

end

%% 2) Correct for wrongly detected gaps
% ensure that every gap is matching a real gap
for i_talker =1:length(CH)

    % check that gap end is associated with an IPU start
    rmIdx=[];
    for i_gap = 1:length(CH{i_talker}.gap.endIdx)
        if isempty(intersect(CH{i_talker}.gap.endIdx(i_gap),CH{i_talker}.IPU.startIdx))
            rmIdx = [rmIdx i_gap];
        end

        %and turn start
        if isempty(intersect(CH{i_talker}.gap.endIdx(i_gap),CH{i_talker}.turn.startIdx))
            rmIdx = [rmIdx i_gap];
        end
  %and turn end
  if isempty(intersect(CH{i_talker}.gap.startIdx(i_gap),CH{CH{i_talker}.gap.channel(i_gap)}.turn.endIdx))
            rmIdx = [rmIdx i_gap];
        end
    end

    CH{i_talker}.gap.startIdx(rmIdx)=[];
    CH{i_talker}.gap.endIdx(rmIdx)=[];
    CH{i_talker}.gap.duration(rmIdx)=[];
    CH{i_talker}.gap.channel(rmIdx)=[];

    CH{i_talker}.gap.numData = numel(CH{i_talker}.gap.duration);

end

%% 3)correct overlapping overlapB and overlapW for same speaker --> multiple
% detections of overlapB for the same speaker --> remove overlapB
for i_talker =1:length(CH)

    % check if start of X and X+1 overlapB are shorter apart then the duration of X
    overlappingWandB = intersect(CH{i_talker}.overlapB.startIdx , CH{i_talker}.overlapW.startIdx);
    overlap_part3{i_talker}=CH{i_talker}.overlapB;
    % if overlapping overlapB's exist
    if ~isempty(overlappingWandB)

        commClassCorrection.multiOverlapBW = length(overlappingWandB);
        for i_overlapWB = 1:length(overlappingWandB)
            %             fprintf('BW_')
            % get overlapB index
            tmp_overlapWBindx=find(CH{i_talker}.overlapB.startIdx == overlappingWandB(i_overlapWB));

            %check if it is a overlapB or overlap W
            if any(CH{i_talker}.overlapB.startIdx(tmp_overlapWBindx)<CH{CH{i_talker}.overlapB.channel(tmp_overlapWBindx)}.turn.endIdx &...
                    CH{i_talker}.overlapB.endIdx(tmp_overlapWBindx)==CH{CH{i_talker}.overlapB.channel(tmp_overlapWBindx)}.turn.endIdx)
                     
                rmIdx=find(CH{i_talker}.overlapW.startIdx==overlappingWandB(i_overlapWB));
                %delete overlapW
                CH{i_talker}.overlapW.startIdx(rmIdx)=[];
            CH{i_talker}.overlapW.endIdx(rmIdx)=[];
            CH{i_talker}.overlapW.duration(rmIdx)=[];
            CH{i_talker}.overlapW.channel(rmIdx)=[];

            CH{i_talker}.overlapW.numData = numel(CH{i_talker}.overlapW.duration);
            else
            % delete overlapping overlap
            CH{i_talker}.overlapB.startIdx(tmp_overlapWBindx)=[];
            CH{i_talker}.overlapB.endIdx(tmp_overlapWBindx)=[];
            CH{i_talker}.overlapB.duration(tmp_overlapWBindx)=[];
            CH{i_talker}.overlapB.channel(tmp_overlapWBindx)=[];

            CH{i_talker}.overlapB.numData = numel(CH{i_talker}.overlapB.duration);
            end

        end
    end

end


%% 4) correct multiple detections of overlapW for the same speaker --> exclude
% on of the two (always identical time and duration)
for i_talker =1:length(CH)
    % check if start of X and X+1 overlapB are shorter apart then the duration of X
    tmp_overlappingDetections = find(diff(CH{i_talker}.overlapW.startIdx) < (CH{i_talker}.overlapW.duration(1:end-1)));

    % if overlapping overlapB's exist
    if ~isempty(tmp_overlappingDetections)
        commClassCorrection.multiOverlapW = length(tmp_overlappingDetections);

        %         fprintf('W')
        % delete overlapping overlap
        CH{i_talker}.overlapW.startIdx(tmp_overlappingDetections+1)=[];
        CH{i_talker}.overlapW.endIdx(tmp_overlappingDetections+1)=[];
        CH{i_talker}.overlapW.duration(tmp_overlappingDetections+1)=[];
        CH{i_talker}.overlapW.channel(tmp_overlappingDetections+1)=[];

        CH{i_talker}.overlapW.numData = numel(CH{i_talker}.overlapW.duration);

    end
end


%% 5) get IPUs after gap and overlapB:
for i_talker =1:length(CH)
    % talkers IPU after gap
    if ~isempty(CH{i_talker}.gap.endIdx)
        [~,idx] = intersect(CH{i_talker}.IPU.startIdx,CH{i_talker}.gap.endIdx);
        CH{i_talker}.IPUafterGap.duration = CH{i_talker}.IPU.duration(idx);
        CH{i_talker}.IPUafterGap.startIdx = CH{i_talker}.IPU.startIdx(idx);
        CH{i_talker}.IPUafterGap.endIdx = CH{i_talker}.IPU.endIdx(idx);
        CH{i_talker}.IPUafterGap.numData = numel(CH{i_talker}.IPU.duration(idx));
    else
        CH{i_talker}.IPUafterGap.duration = [];
        CH{i_talker}.IPUafterGap.startIdx = [];
        CH{i_talker}.IPUafterGap.endIdx = [];
        CH{i_talker}.IPUafterGap.numData = [];
    end

    % talkers IPU after overlap[B]
    if ~isempty(CH{i_talker}.overlapB.endIdx)
        [~,idx] = intersect(CH{i_talker}.IPU.startIdx,CH{i_talker}.overlapB.startIdx);%

        CH{i_talker}.IPUafterOverlap.duration = CH{i_talker}.IPU.duration(idx);
        CH{i_talker}.IPUafterOverlap.startIdx = CH{i_talker}.IPU.startIdx(idx);
        CH{i_talker}.IPUafterOverlap.endIdx = CH{i_talker}.IPU.endIdx(idx);
        CH{i_talker}.IPUafterOverlap.numData = numel(CH{i_talker}.IPU.duration);
    else
        CH{i_talker}.IPUafterOverlap.duration = [];
        CH{i_talker}.IPUafterOverlap.startIdx = [];
        CH{i_talker}.IPUafterOverlap.endIdx = [];
        CH{i_talker}.IPUafterOverlap.numData = [];
    end
end


%% 6) get IPUs before gap and overlapB:
for i_talker =1:length(CH)

    if ~isempty(CH{i_talker}.gap.startIdx)
        for i_gap = 1:CH{i_talker}.gap.numData

            [~,idx] = intersect(CH{CH{i_talker}.gap.channel(i_gap)}.IPU.endIdx, CH{i_talker}.gap.startIdx(i_gap)); %ORG

            CH{i_talker}.IPUbeforeGap.duration(i_gap) = CH{CH{i_talker}.gap.channel(i_gap)}.IPU.duration(idx);
            CH{i_talker}.IPUbeforeGap.startIdx(i_gap) = CH{CH{i_talker}.gap.channel(i_gap)}.IPU.startIdx(idx);
            CH{i_talker}.IPUbeforeGap.endIdx(i_gap) = CH{CH{i_talker}.gap.channel(i_gap)}.IPU.endIdx(idx);
        end

        if CH{i_talker}.IPUbeforeGap.numData > i_gap
            CH{i_talker}.IPUbeforeGap.duration(i_gap+1:end) = [];
            CH{i_talker}.IPUbeforeGap.startIdx(i_gap+1:end) = [];
            CH{i_talker}.IPUbeforeGap.endIdx(i_gap+1:end) = [];
        end
    else
        CH{i_talker}.IPUbeforeGap.duration = [];
        CH{i_talker}.IPUbeforeGap.startIdx = [];
        CH{i_talker}.IPUbeforeGap.endIdx= [];
    end
    CH{i_talker}.IPUbeforeGap.numData = numel(CH{i_talker}.IPUbeforeGap.endIdx);


    if ~isempty(isempty(CH{i_talker}.overlapB.startIdx))
        if CH{i_talker}.overlapB.numData>0 % in the case that no overlap is found --> set IPUbefore to zero
            for i_ovB = 1:CH{i_talker}.overlapB.numData
                [~,idx] = intersect(CH{CH{i_talker}.overlapB.channel(i_ovB)}.IPU.endIdx,  CH{i_talker}.overlapB.startIdx(i_ovB) : CH{i_talker}.overlapB.endIdx(i_ovB));

                CH{i_talker}.IPUbeforeOverlap.duration(i_ovB) = CH{CH{i_talker}.overlapB.channel(i_ovB)}.IPU.duration(idx);
                CH{i_talker}.IPUbeforeOverlap.startIdx(i_ovB) = CH{CH{i_talker}.overlapB.channel(i_ovB)}.IPU.startIdx(idx);
                CH{i_talker}.IPUbeforeOverlap.endIdx(i_ovB) = CH{CH{i_talker}.overlapB.channel(i_ovB)}.IPU.endIdx(idx);

            end


        else
            i_ovB=0;
        end
        % make sure that old values are removed
        if CH{i_talker}.IPUbeforeOverlap.numData > i_ovB
            CH{i_talker}.IPUbeforeOverlap.duration(i_ovB+1:end) = [];
            CH{i_talker}.IPUbeforeOverlap.startIdx(i_ovB+1:end) = [];
            CH{i_talker}.IPUbeforeOverlap.endIdx(i_ovB+1:end) = [];
        end
    else
        CH{i_talker}.IPUbeforeOverlap.duration = [];
        CH{i_talker}.IPUbeforeOverlap.startIdx = [];
        CH{i_talker}.IPUbeforeOverlap.endIdx= [];
    end

    CH{i_talker}.IPUbeforeOverlap.numData = numel(CH{i_talker}.IPUbeforeOverlap.endIdx);
end

%% Detect first channel
     idx=find(arrayfun(@(x) CH{x}.turn.numData>0,1:length(CH)));
        [~,firstTurn]=min(arrayfun(@(x) CH{x}.turn.startIdx(1),idx));
        firstTurn=idx(firstTurn);
            CH{firstTurn}.overlapB.startIdx=[1,CH{firstTurn}.overlapB.startIdx];
            CH{firstTurn}.overlapB.endIdx=[1,CH{firstTurn}.overlapB.endIdx];
            CH{firstTurn}.overlapB.duration=[NaN,  CH{firstTurn}.overlapB.duration];
            CH{firstTurn}.overlapB.channel=[NaN,CH{firstTurn}.overlapB.channel];
            
            
%% SANITY CHECK...

% if no of iPUs after gap does not match
for i_talker =1:length(CH)
    if ne(length(CH{i_talker}.gap.startIdx),length(CH{i_talker}.IPUafterGap.startIdx))
        fprintf('SOMETHING IS WRONG  GAPafter')
        keyboard
    end


    % if no of iPUs before gap does not match
    if ne(length(CH{i_talker}.gap.startIdx),length(CH{i_talker}.IPUbeforeGap.startIdx))
        fprintf('SOMETHING IS WRONG GAPbefore ')
        keyboard
    end

    % if no of iPUs after overlap does not match
  

    if ne(CH{i_talker}.overlapB.numData,length(CH{i_talker}.IPUafterOverlap.startIdx))
        fprintf('SOMETHING IS WRONG overlapafter')
        keyboard
    end

    if ne(CH{i_talker}.overlapB.numData,length(CH{i_talker}.IPUbeforeOverlap.startIdx))
        fprintf('SOMETHING IS WRONG overlapbefore')
        keyboard
    end

    if ne(CH{i_talker}.overlapB.numData+CH{i_talker}.gap.numData,CH{i_talker}.turn.numData)
        %check if it is the first turn
       
        if i_talker ~= firstTurn ||...
                (i_talker == firstTurn && ne(CH{i_talker}.overlapB.numData+CH{i_talker}.gap.numData+1,CH{i_talker}.turn.numData))
        
%             tmpIdx=find(~ismember(CH{i_talker}.turn.startIdx,[CH{i_talker}.overlapB.startIdx,CH{i_talker}.gap.endIdx]));
% potential_gaps=arrayfun(@(x) find(CH{i_talker}.turn.startIdx(tmpIdx)>=CH{x}.turn.endIdx&[CH{i_talker}.turn.startIdx(tmpIdx)<CH{x}.turn.startIdx(2:end),1]),idx,'UniformOutput',false);
% tmpchannel=[];
% tmpgapStart=[];
% for i=idx
%     if ~isempty(potential_gaps{i})
%         tmpchannel(ii)=repmat(i,size(potential_gaps{i}));
%         tmpgapStart=CH{i}.turn.endIdx(potential_gaps{i})
%     end
% end
% if ~isempty(tmpchannel)
%     [gapStart,maxIdx]=max(tmpgapStart);
%     CH{i_talker}.gap.startIdx=[CH{i_talker}.gap.startIdx,gapStart];
%     CH{i_talker}.gap.endIdx=[CH{i_talker}.gap.endIdx,CH{i_talker}.turn.startIdx(tmpIdx)];
%     CH{i_talker}.gap.duration=[CH{i_talker}.gap.duration,CH{i_talker}.turn.startIdx(tmpIdx)-gapStart];
%     CH{i_talker}.gap.channel=[CH{i_talker}.gap.channel,tmpchannel(maxIdx)];
% 
%     CH{i_talker}.gap.numData = numel(CH{i_talker}.gap.duration);
% %check if gap
%             if CH{i_talker}.turn.startIdx(tmpIdx)
                    fprintf('SOMETHING IS WRONG number FTOs')

          % keyboard
        end
    end
end


CH_reorg = CH;
end