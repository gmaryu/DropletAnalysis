function oscillation_table=addTuningInt(oscillation_table, sumTrackTable, aquisitionPosition, DyeColumn)
%% add tuning information on mlx file
%{
Description of function
    Add dye intensity at 1st frame to oscillation_table

oscillation_table: all oscillation information
sumTrackTable: droplet tracking results
aquisitionPosition: target positions
DyeColumn: the name of column of dye channel
%}

clear positions 
clear trackIndex
clear tracks
clear track
clear oscis
clear dropletID
clear positionIndex
%aquisitionPosition=genPos;
%DyeColumn='DAPI_whole';

positions=extractfield(oscillation_table,'position');
if length(positions)~=length(sumTrackTable)
    disp('Length of osic_table and sumTrackTable is different.');
    return
end

positionIndex=1;
while positionIndex <= length(positions)
    % tracking data
    trackIndex=aquisitionPosition==positions(positionIndex);
    tracks=sumTrackTable{trackIndex};
    dropletIds=extractfield(tracks, 'id');
    
    % oscillation data
    oscis=oscillation_table(positionIndex).data;
    
    % main
    for o = 1:size(oscis,1)
        dropletID=oscis.dropID(o);
        dropletIndex=dropletIds==dropletID;
        track=tracks(dropletIndex).Feat;
        eval(['tuningDyeIntensity=track.',DyeColumn,'(1);']);
        oscis.DyeIntensity(o)=tuningDyeIntensity;
    end
       
    %update data
    oscillation_table(positionIndex).data=oscis;
    positionIndex=positionIndex+1;
end
end
