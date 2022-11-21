function f=dropletPeriods_box(oscillation_table, pos)
%{
Description
    plor a boxplot for designated positions
Args
    oscillation_table
    pos
%}
%%
%osci_table=oscillation_table;
%pos=genPos;
f=figure();
hold on
fidx=[];
pnum=[];
x_box=[];
g_box=[];
for n=1:length(oscillation_table)
   fidx=[fidx n]; % field index 
   pnum=[pnum oscillation_table(n).position];
end

for i=1:length(pos)
    f=find(pnum==pos(i));
    PosN_lastpeaks=[];
    PosN_osci=osci_table(fidx(f)).data;
    PosN_unq_dropletID=unique(PosN_osci.dropID);
    periodpeak=PosN_osci.periodPeak;
    
    x_box=[x_box; periodpeak]; % last peak timepoint data
    tmp_g=repmat({['Pos',int2str(pos(i))]},length(periodpeak),1); % make a group index for plotting
    g_box=[g_box; tmp_g];
    
end
boxplot(x_box,g_box);
ylabel('Period(Frame number)');
%hold off
%savefig([savePath,'/boxPlot_LastPeaks_test.fig']);
end