function sumTrackTable=skewFindChangepts(sumTrackTable)
%%
for p=1:size(sumTrackTable,2)
    tracks=sumTrackTable{p};
    for i=1:size(tracks)
        tmpSkewness=tracks(i).Feat.skewness;
        findchangepts(tmpSkewness,'MaxNumChanges',8);
        fname = ['TrackID_',num2str(tracks(i).id),'.png'];
        saveas(gcf, fname);
        %ipt=findchangepts(tmpSkewness,'MaxNumChanges',8);
    end    
end

