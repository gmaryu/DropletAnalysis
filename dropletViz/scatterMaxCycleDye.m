function scatterMaxCycleDye(oscillation_table, position)

allPosition=extractfield(oscillation_table, 'position');
selected_position=ismember(allPosition, position);

PosX_maxcycleID=[];
PosX_dye=[];

for i=1:length(selected_position)
    if selected_position(i)==1
        tmp_table=oscillation_table(i).data;
        col_dropletID=tmp_table.dropID;
        tmp_unique_d=unique(col_dropletID);
        col_dye=tmp_table.DyeIntensity;
        col_cycle=tmp_table.cycleID;
        
        for j=1:length(tmp_unique_d)
            % cycle info
            did=tmp_unique_d(j);
            
            % droplet dye info
            dropInd=find(col_dropletID==did);
            %disp(dropInd);
            dyeInt=max(col_dye(dropInd));
            max_cycnum=max(col_cycle(dropInd));
            
            PosX_maxcycleID=[PosX_maxcycleID, max_cycnum];
            PosX_dye=[PosX_dye,dyeInt];
            
        end
    end
end

scatter(PosX_dye, PosX_maxcycleID,15,'filled','MarkerFaceColor','#646464');
xlabel('Dye Intensity');
ylabel('Cycle Number');

end