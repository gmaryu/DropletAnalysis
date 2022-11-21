function [tbl_all, tbl_CN]=swarmPlot_npeaks(oscillation_table, positions, option, ratio_thresh)
%{
Description
    plot a bar plot for the number of peaks (including 1st peak) in each
    position
Args
    pos: position what you plot
    osci_table: oscillation data table
    group: group of experiment condition. assign as integer
        e.g. [1,1,1,2,2,2,3,3,3,4,4,4,4,4]
%}

%% 
% 
% positions=[6,7,8];
% option=2;
% ratio_thresh=0.1;

%%

tbl_all=table;
tbl_CN=table;

tlabel={};
tlabel_CN={};


allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, positions);

for s=1:length(selected_position)
    if selected_position(s)==1
        PosN_osci=oscillation_table(s).data;
        PosN_unq_dropletID=unique(PosN_osci.dropID);
        PosN={};
        npeaks_data=[];
        PosN_nuc={};
        PosN_cyto={};
        periodpeak=PosN_osci.periodPeak;
        cycleID=PosN_osci.cycleID;
        nuc_seg_ratio=PosN_osci.seg_ratio;
        dropletID=PosN_osci.dropID;
        
        for i=1:length(PosN_unq_dropletID)
            PosN=vertcat(PosN,['Pos ',int2str(allpositions(s))]);
            npeaks=length(find(PosN_osci.dropID==PosN_unq_dropletID(i))); % number of peaks
            npeaks_data=[npeaks_data, npeaks];
            tmp_nuc_seg_ratio=nuc_seg_ratio(PosN_osci.dropID==PosN_unq_dropletID(i));
            if max(tmp_nuc_seg_ratio) > ratio_thresh
                PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' N']);
                tmp_nuc={['Pos ',int2str(allpositions(s)),' N'],npeaks,PosN_unq_dropletID(i)};
                tbl_CN=vertcat(tbl_CN,tmp_nuc);
                tlabel_CN=horzcat(tlabel_CN, ['Pos ',int2str(allpositions(s)),' N']);
            else
                PosN_nuc=vertcat(PosN,['Pos ',int2str(allpositions(s)),' C']);
                tmp_cyto={['Pos ',int2str(allpositions(s)),' C'],npeaks,PosN_unq_dropletID(i)};
                tbl_CN=vertcat(tbl_CN,tmp_cyto);
                tlabel_CN=horzcat(tlabel_CN, ['Pos ',int2str(allpositions(s)),' C']);
            end
        end        
            
        tmp_tbl=table(PosN,npeaks_data');
        tbl_all=vertcat(tbl_all,tmp_tbl);
        tlabel=horzcat(tlabel, ['Pos ',int2str(allpositions(s))]);
        x=categorical(tbl_all.PosN, tlabel);
        y = tbl_all.Var2;
        c= tbl_all.Var2;
        
        
        if exist('tbl_CN') && ~isempty(tbl_CN) 
            tbl_CN.Properties.VariableNames{1} = 'PosN';
            tbl_CN.Properties.VariableNames{2} = 'MaxCycleID';
            tbl_CN.Properties.VariableNames{3} = 'dropletID';
        end
        

    end
end

figure();
swarmchart(x,y,20,c,'filled');
ylabel('Period (min)')

cbr=colorbar;
cbr.Label.String='Cycle Number';

if option==2
    figure();
    x1=categorical(tbl_CN.PosN);
    y1=tbl_CN.MaxCycleID;
    c1=tbl_CN.MaxCycleID;
    figure();
    swarmchart(x1,y1,20,c1,'filled');
    ylabel('Period (min)')

    cbr=colorbar;
    cbr.Label.String='Cycle Number';
end
