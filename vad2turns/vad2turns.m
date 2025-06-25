function [mVADout,t,vspstate,tovl] = vad2turns( vT, mVAD, maxgapdur, minsegmentdur )
    if nargin < 3
        maxgapdur = 0.18;
        % value form Heldner2010
    end
    if nargin < 4
        minsegmentdur = 0.09;
        % value from Heldner2010
    end
    mVADout = zeros(size(mVAD));
    cArg = {};
    for ch=1:size(mVAD,2)
        act = mVAD(:,ch);
        act(1) = 0;
        act(end+1) = 0;
        dact = diff(act);
        tstart = vT(find(dact==1));
        tend = vT(find(dact==-1));
        [tstart, tend] = remove_gaps( tstart, tend, maxgapdur );
        % length filtering:
        [tstart,tend] = remove_short_segments( tstart, tend, minsegmentdur );
        for k=1:numel(tstart)
            idx_start = find(vT==tstart(k),1);
            idx_end = find(vT==tend(k),1);
            mVADout(idx_start:idx_end,ch) = 1;
        end
        cArg{end+1} = tstart;
        cArg{end+1} = tend;
    end
    mStateChange = diff(mVADout,1,1);
    [t,vspstate,tovl] = get_turntakes( cArg{:} );
end
