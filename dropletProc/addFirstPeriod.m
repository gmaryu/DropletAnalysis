function oscillation_table=addFirstPeriod(oscillation_table, sumTrackTable)
% add 1st period
%%
clear positions 
clear trackIndex
clear oscis
clear dropletID
clear positionIndex
%aquisitionPosition=genPos;
%DyeColumn='DAPI_whole';

positions=extractfield(oscillation_table,'position');
if length(positions)~=length(sumTrackTable)
    disp('Length of osic_table and sumTrackTable is different.');
end

positionIndex=1;
while positionIndex <= length(positions)    
    % oscillation data
    oscis=oscillation_table(positionIndex).data;
    
    % main
    for o = 1:size(oscis,1)
        dropletID=oscis.dropID(o);
        samedroplets=oscis.dropID==dropletID;
        firstPeriod=nanmin(oscis.periodPeak(samedroplets));
        oscis.firstPeriod(o)=firstPeriod;
    end
       
    % update data
    oscillation_table(positionIndex).data=oscis;
    positionIndex=positionIndex+1;
end
end