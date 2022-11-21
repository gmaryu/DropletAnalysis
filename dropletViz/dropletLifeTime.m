function f=dropletLifeTime(osci_table, pos, savePath)
%{
Description
    plot a box plot for designated positions
Args
    osci_table: oscillation table
    pos: a list of positions
    savePath: data save path
%}

f=figure();
hold on
fidx=[];
pnum=[];
x_box=[];
g_box=[];
for n=1:length(osci_table)
   fidx=[fidx n]; % field index 
   pnum=[pnum osci_table(n).position];
end

for i=1:length(pos)
    f=find(pnum==pos(i));
    PosN_lastpeaks=[];
    PosN_osci=osci_table(fidx(f)).data;
    PosN_unq_dropletID=unique(PosN_osci.dropID);
    peaktime=PosN_osci.peakTime;
    
    for j=1:length(PosN_unq_dropletID)
        id=PosN_unq_dropletID(j);
        tmpdroplet=find(PosN_osci.dropID==id);
        PosN_lastpeaks(j)=peaktime(tmpdroplet(end));
    end
    x_box=[x_box; PosN_lastpeaks']; % last peak timepoint data
    tmp_g=repmat({['Pos',int2str(pos(i))]},j,1); % make a group index for plotting
    g_box=[g_box; tmp_g];
    
end
boxplot(x_box,g_box);
ylabel('Frame number');
hold off
savefig([savePath,'/boxPlot_LastPeaks_test.fig']);
end