function oscillation_table=dyeint2molar2(sumTrackTable, aquisitionPosition,DyeColumn,max_M,min_M)
% % % %
% Description of function
%    convert intensity data to molar data
% 
% sumTrackTable: generate struture by tracking result
% aquisitionPosition: list of position
% DyeColumn: colum name of dye intensity of sumTrackTable
% max_M: maximum molar concentration
% min_M: minumum molar concentration
% % % %

clear positions 
clear trackIndex
clear tracks
clear track
clear oscis
clear dropletID
clear positionIndex

aquisitionPosition=genPos;
DyeColumn='DAPI_whole';
max_M=10;
min_M=0;

% add dye intensity to oscillation table
positionIndex=1;
while positionIndex <= length(positions)
    % tracking data
    trackIndex=aquisitionPosition==positions(positionIndex);
    tracks=sumTrackTable{trackIndex};
    
    % oscillation data
    oscis=oscillation_table(positionIndex).data;
    
    % main
    for o = 1:size(oscis,1)
        dropletID=oscis.dropID(o);
        track=tracks(dropletID).Feat;
        eval(['tuningDyeIntensity=track.',DyeColumn,'(1);']);
        oscis.DyeIntensity(o)=tuningDyeIntensity;
    end
       
    %update data
    oscillation_table(positionIndex).data=oscis;
    positionIndex=positionIndex+1;
end

% convert intensity data to molar information
% extract max, min intensity
positionIndex=1;
droplets_int=[];
while positionIndex <= length(positions)
    % oscillation data
    oscis=oscillation_table(positionIndex).data;
    % 
    for o = 1:size(oscis,1)
        tmp_dye=oscis.DyeIntensity;
        droplets_int=[droplets_int;tmp_dye];
    end
    positionIndex=positionIndex+1;
end
droplets_int=unique(droplets_int);
upper_int=prctile(droplets_int,98);
lower_int=prctile(droplets_int,2);

%%
%disp(max_int)
% calculate estimated molar value
positionIndex=1;
while positionIndex <= length(positions)
    % tracking data
    trackIndex=aquisitionPosition==positions(positionIndex);
    tracks=sumTrackTable{trackIndex};
    
    % oscillation data
    oscis=oscillation_table(positionIndex).data;
    oscis.est_Molar=(oscis.DyeIntensity-lower_int)/(upper_int-lower_int)*(max_M-min_M);
    
    oscis.est_Molar(oscis.est_Molar>max_M)=max_M;
    oscis.est_Molar(oscis.est_Molar<min_M)=min_M;
    
    %update data
    oscillation_table(positionIndex).data=oscis;
    positionIndex=positionIndex+1;
end

end

