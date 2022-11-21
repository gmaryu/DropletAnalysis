function tracks=pbCorrection(tracks,colNames)
%{
tracks: tracking result object for target position
colNames: list of column names 
        {1}:reference column e.g. *_ratio_whole
        {2}:target column e.g. *_ratio_nuc
%}
    %% mean value of target channel
    tmp=NaN(size(tracks,1),max([tracks.Len].')); %array for channel of interest
    for n =1:size(tracks, 1)
        eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},';']);
        tmp(n,1:length(tmpSig))=tmpSig;
    end
    mtmp=nanmean(tmp,1); % vector of mean intensity of channel of interest
    
    %% fitting function and estimated value
    [p,s,mu] = polyfit((1:numel(mtmp)),mtmp,6); % fitting
    f_y = polyval(p,(1:numel(mtmp)),[],mu); % estimated vector
    
    %% data correction by estimated value
    if length(colNames)==1    
        for n =1:size(tracks)
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},';']);
            ntmp=tmpSig./f_y(1:length(tmpSig))';
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_pb=ntmp;']);
        end
    else
        for n =1:size(tracks)
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{2},';']);
            ntmp=tmpSig./f_y(1:length(tmpSig))';
            eval(['tracks(',num2str(n),').Feat.',colNames{2},'_pb=ntmp;']);
        end
    end
end
