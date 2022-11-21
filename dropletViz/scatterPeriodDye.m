function scatterPeriodDye(oscillation_table, sumTrackTable, position)
%{
    oscillation_table: dye intensity should be added
    sumTrackTable: 
    position: position list to be plotted
%}

%%
allPosition=extractfield(oscillation_table, 'position');
selected_position=ismember(allPosition, position);

%%
PosX_period=[];
PosX_cycleID=[];
PosX_dye=[];

for i=1:length(selected_position)
    if selected_position(i)==1
        %disp(allPosition(i));
        PosN_osci=oscillation_table(i).data;
        PosN_tracks=sumTrackTable{i};
        PosN_tracks_did=extractfield(PosN_tracks,'id');
        dropletID=PosN_osci.dropID;
        periodpeak=PosN_osci.periodPeak;
        peaktime=PosN_osci.peakTime;
        cycleID=PosN_osci.cycleID;
        
        for j=1:length(dropletID)
            % cycle info
            did=dropletID(j);
            cid=cycleID(j);
            pp=periodpeak(j);
            pt=peaktime(j);
            % droplet dye info
            dropInd=find(PosN_tracks_did==did);
            %disp(dropInd);
            dyeInt=PosN_tracks(dropInd).Feat.BFP_whole(1);
            
            % make vector for plot
            PosX_cycleID=vertcat(PosX_cycleID, cid);
            PosX_dye=vertcat(PosX_dye, dyeInt);
            if isnan(pp) && cid==1
                PosX_period=vertcat(PosX_period, pt);
            else
                PosX_period=vertcat(PosX_period, pp);
            end
            
        end
    
    end    
    
end

scatter(PosX_dye, PosX_period,5,PosX_cycleID,'filled');
colorbar;
xlabel('Dye Intensity');
ylabel('Period length (Frame)');

end