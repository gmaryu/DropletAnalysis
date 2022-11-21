function plotDyeVSCycleNumber(position, oscillation_table, sumTrackTable)
%{
    plot a swarmplot for "Tuned dye intensity" vs "cycle number(inclue one peak)"
    position: list of position where you want to analyze
    oscillation_table: a table of oscillation information
    sumTrackTable: a structure of summary of tracking results
%}

%% 
allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, position);

if length(sumTrackTable)~=length(oscillation_table)
    disp('Number of position of Tracking results and  is mismatch.');
end

%% get max cycle number
max_cycle_number=0;
for p =1:length(oscillation_table)
    if selected_position(p)==1
        tmp_table=oscillation_table(p).data;
        tmp_cycleIDs=tmp_table.cycleID;
        %disp(max(tmp_cycleIDs));
        if max(tmp_cycleIDs) > max_cycle_number
            max_cycle_number = max(tmp_cycleIDs);
        end
    end
end
disp(['Max cycle number is ', num2str(max_cycle_number)]);

%% generate variables
DyeInt0=[]; % for no ocsillation droplets
for n = 1:max_cycle_number
    %disp(n);
    eval(['DyeInt',num2str(n),'=[];']);
end

%% add dye intensity values from oscillation_table
for p =1:length(oscillation_table)
    if selected_position(p)==1
        tmp_position=oscillation_table(p).position;
        tmp_table=oscillation_table(p).data;
        tmp_cycleIDs=tmp_table.cycleID;
        tmp_dyeIntensity=tmp_table.DyeIntensity;
        tmp_dropletID=tmp_table.dropID;
        tmp_unq_dropletID=unique(tmp_table.dropID);
        
        for id = 1:length(tmp_unq_dropletID)
            tmp_drop=tmp_unq_dropletID(id);
            tmp_dye=mean(tmp_dyeIntensity(tmp_dropletID==tmp_drop));
            tmp_max_cycleNum=max(tmp_cycleIDs(tmp_dropletID==tmp_drop));
            eval(['DyeInt',num2str(tmp_max_cycleNum),'=[DyeInt',num2str(tmp_max_cycleNum),', tmp_dye];']);
        end
    end
end

%% add dye intensity values of no-oscillation droplets
for p =1:length(oscillation_table)
    if selected_position(p)==1
        
        tmp_tracks=sumTrackTable{p};
        tmp_all_tracked_id=extractfield(tmp_tracks, 'id');
        
        tmp_position=oscillation_table(p).position;
        tmp_table=oscillation_table(p).data;
        tmp_dropletID=tmp_table.dropID;
        tmp_unq_dropletID=unique(tmp_table.dropID);
        
        for id = 1:length(tmp_all_tracked_id)
            if ~ismember(tmp_unq_dropletID, tmp_all_tracked_id(id))
                tmp_droplet=tmp_tracks(id).Feat;
                tmp_v = tmp_droplet.BFP_whole(1);
                DyeInt0=[DyeInt0, tmp_v]; 
                %disp('no oscillation')
            end
        end
    end
end


%% plot figures
clear gca;
figure();
x0=0*ones(1,length(DyeInt0));
swarmchart(x0, DyeInt0,20,'filled');
hold on
for n = 1:max_cycle_number
    if eval(['isempty(DyeInt',num2str(n),')'])
        eval(['DyeInt',num2str(n),'=NaN;'])
    end
    eval(['x',num2str(n),'=',num2str(n),'*ones(1,length(DyeInt',num2str(n),'));'])
    eval(['swarmchart(x',num2str(n),', DyeInt',num2str(n),',20, "filled");'])
end
ylabel('Dye Intensity');
xlabel('Number of peaks')
xlim([-0.5,max_cycle_number+0.5]);
hold off
