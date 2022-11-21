function plotPeriod_swarmChart2(oscillation_table, positions, groups, interval_time, option, cycleid, yrange)
%{
Description

Args
    oscillation_table: 
    positions: positions of interest
    groups: a list of condition
        e.g.: [1,1,1,2,2,2,3,3,3,3,3,3]
    interval_time: the value of interval time (min)
    yrange: a list for ylim function
    options: 
        1 -> all
        2 -> all, nucleus, cytoplasm
    cycleid: cycle id of interest
        e.g: [2,3,4] (first 3 whole oscillation periods)
%}
%%
% positions = [3,4,5,6,7,8,9,10,11,12,13,14,15,16];
% groups = [1,1,1,2,2,2,3,3,3,3,3,3,3,3];
% interval_time = 4;
% option = 2;
% cycleid = [2,3,4];

clear gca
%%

tbl_all=table;
tbl_CN=table;

tlabel={};
tlabel_CN={};

allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, positions);


for s=1:length(selected_position)
    %disp(s)
    if selected_position(s)==1
        tmp_position=allpositions(s); % position of interest
        tmp_group=groups(ismember(positions,allpositions(s))); % group tag of position of interest
        PosN_osci=oscillation_table(s).data; % data chunk of designated position
        PosN_unq_dropletID=unique(PosN_osci.dropID);
        PosN={};
        PosN_nuc={};
        PosN_cyto={};
        tmp_cycleID=[];
        tmp_period=[];
        tmp_groups=[];
        
        periodpeak=PosN_osci.periodPeak;
        cycleID=PosN_osci.cycleID;
        nuc_seg=PosN_osci.Nuc;
        dropletID=PosN_osci.dropID;
        
        for i=1:length(periodpeak)
            if max(ismember(cycleid, cycleID(i)))==1
                PosN=vertcat(PosN,['Pos ',int2str(allpositions(s))]);
                if nuc_seg(i)==1
                    PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' N']);
                    tmp_nuc={['Pos ',int2str(allpositions(s)),' N'], periodpeak(i), cycleID(i), dropletID(i), ['Group ', int2str(tmp_group),' N']};
                    tbl_CN=vertcat(tbl_CN,tmp_nuc);
                    tlabel_CN=horzcat(tlabel_CN, ['Group ', int2str(tmp_group),' N']);
                elseif nuc_seg(i)==0
                    PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' C']);
                    tmp_cyto={['Pos ',int2str(allpositions(s)),' C'], periodpeak(i), cycleID(i), dropletID(i), ['Group ', int2str(tmp_group),' C']};
                    tbl_CN=vertcat(tbl_CN,tmp_cyto);
                    tlabel_CN=horzcat(tlabel_CN, ['Group ', int2str(tmp_group),' C']);
                else
                    disp('!!!')
                end
                tmp_period=[tmp_period; periodpeak(i)];
                tmp_cycleID=[tmp_cycleID; cycleID(i)];
                tmp_groups=[tmp_groups;['Group ', int2str(tmp_group)]];
            end
        end
        tmp_tbl=table(PosN,tmp_period,tmp_cycleID,cellstr(tmp_groups));
        tbl_all=vertcat(tbl_all,tmp_tbl);
        tlabel=horzcat(tlabel, ['Group ',int2str(tmp_group)]);
        x=categorical(tbl_all.Var4);
        y = tbl_all.tmp_period*interval_time;
        c= tbl_all.tmp_cycleID;
        
        
        if exist('tbl_CN') && ~isempty(tbl_CN) 
            tbl_CN.Properties.VariableNames{1} = 'PosN';
            tbl_CN.Properties.VariableNames{2} = 'periodpeak';
            tbl_CN.Properties.VariableNames{3} = 'cycleID';
            tbl_CN.Properties.VariableNames{4} = 'dropletID';
            tbl_CN.Properties.VariableNames{5} = 'group';
        end
        

    end
end
%%
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
    x1=categorical(tbl_CN.group);
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