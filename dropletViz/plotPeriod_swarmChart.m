function [tbl_all, tbl_CN]=plotPeriod_swarmChart(oscillation_table, positions, interval_time, yrange, option)
%{
% interval_time: the value of interval time (min)
% yrange: a list for ylim function
% options: 
    1 -> all
    2 -> all, nucleus, cytoplasm
%}
clear gca
%%
posIdx=1;
period_data=[];
period_stdev=[];
period_nuc_data=[];
period_nuc_stdev=[];
period_cyto_data=[];
period_cyto_stdev=[];

tbl_all=table;
tbl_CN=table;

tlabel={};
tlabel_CN={};

fidx=[];
pnum=[];

allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, positions);

for s=1:length(selected_position)
    if selected_position(s)==1
        PosN_osci=oscillation_table(s).data;
        PosN_unq_dropletID=unique(PosN_osci.dropID);
        PosN={};
        PosN_nuc={};
        PosN_cyto={};
        periodpeak=PosN_osci.periodPeak;
        cycleID=PosN_osci.cycleID;
        nuc_seg=PosN_osci.Nuc;
        dropletID=PosN_osci.dropID;
        for i=1:length(periodpeak)
            PosN=vertcat(PosN,['Pos ',int2str(allpositions(s))]);
            if nuc_seg(i)==1
                PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' N']);
                tmp_nuc={['Pos ',int2str(allpositions(s)),' N'],periodpeak(i),cycleID(i),dropletID(i)};
                tbl_CN=vertcat(tbl_CN,tmp_nuc);
                tlabel_CN=horzcat(tlabel_CN, ['Pos ',int2str(allpositions(s)),' N']);
            else
                PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' C']);
                tmp_cyto={['Pos ',int2str(allpositions(s)),' C'],periodpeak(i),cycleID(i),dropletID(i)};
                tbl_CN=vertcat(tbl_CN,tmp_cyto);
                tlabel_CN=horzcat(tlabel_CN, ['Pos ',int2str(allpositions(s)),' C']);
            end
        end
        tmp_tbl=table(PosN,periodpeak,cycleID);
        tbl_all=vertcat(tbl_all,tmp_tbl);
        tlabel=horzcat(tlabel, ['Pos ',int2str(allpositions(s))]);
        x=categorical(tbl_all.PosN, tlabel);
        y = tbl_all.periodpeak*interval_time;
        c= tbl_all.cycleID;
        
        
        if exist('tbl_CN') && ~isempty(tbl_CN) 
            tbl_CN.Properties.VariableNames{1} = 'PosN';
            tbl_CN.Properties.VariableNames{2} = 'periodpeak';
            tbl_CN.Properties.VariableNames{3} = 'cycleID';
            tbl_CN.Properties.VariableNames{4} = 'dropletID';
        end
        

    end
end

figure();
swarmchart(x,y,20,c,'filled');
ylabel('Period (min)')
if ~isempty(yrange)
    ylim(yrange)
end
cbr=colorbar;
cbr.Label.String='Cycle Number';

if option==2
    figure();
    x1=categorical(tbl_CN.PosN);
    y1=tbl_CN.periodpeak*interval_time;
    c1=tbl_CN.cycleID;
    figure();
    swarmchart(x1,y1,20,c1,'filled');
    ylabel('Period (min)')
    if ~isempty(yrange)
        ylim(yrange)
    end
    cbr=colorbar;
    cbr.Label.String='Cycle Number';
end

%savefig('barPlot_npeaks.fig')
%end