function dropletXY(sumTrackTable, genPos)
%{
Description
    Plot a moving trajectory of indivisual droplet
Args
    sumTrackTable: tracking result matlab structure
    genPos: Position of interest, 'genPos' is for all positions
%}
%%
for i=1:length(genPos)
    tracks=sumTrackTable{i};
    x_pos=[];
    y_pos=[];
    cmp=[];
    %cmap=jet(numel([1:length(tracks)])');
    for j=1:length(tracks)
       cmp(j)=j;
       x_pos(j)=tracks(j).Feat.xcoord(1);
       y_pos(j)=tracks(j).Feat.ycoord(1);
    end
    figure();
    pointsize=10;
    scatter(x_pos, y_pos,pointsize,cmp,'filled');
    colorbar();
    fn=['xyPlot_1stFrame_Pos',int2str(genPos(i)),'.fig'];
    savefig(fn);
end
end