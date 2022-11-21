function oscillation_table=addNucInfo2OsciTable(oscillation_table, sumTrackTable, qposition,thresh)
%{
Description
    add nuclear segmentation information to oscillation_table
    1. data.Nuc: nucleus was detected in a certain droplet (logical)
    2. data.duration: quantification initial point to first peak and peak to peak
                frame number after 2nd peaks
    3. data.seg_frame: nucleus detected frame number in a certain
    oscillation (integer)
    3. data.seg_ratio: ratio of nucleus detected ratio in a certain
    oscillation. "data.seg_frame"./"data.period" (double)
  
Args 
    oscillation_table: genarated table by peak detection function
    sumTrackTable: generate struture by tracking result
    aquisitionPosition: list of positions to analyze
    thresh: percentage of nuclear detected frame in a certain droplet (int)
%}
%%
clear positions 
clear trackIndex
clear tracks
clear track
clear oscis
clear dropletID
clear positionIndex
%aquisitionPosition=pos;
%tmpsumTrackTable=sumTrackTable(1:5);
allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, qposition);

if length(allpositions)~=length(sumTrackTable)
    disp('Length of osicllation_table and sumTrackTable is different.');
    return
end
%%



% add dye intensity to oscillation table

for positionIndex = 1:length(allpositions)
    
    if selected_position(positionIndex)==1
        % tracking data
        tracks=sumTrackTable{positionIndex};
        tracks_ids=extractfield(tracks, 'id');
        
        % oscillation data
        oscillations=oscillation_table(positionIndex).data;
        % remove dropID==NaN from oscillation_table
        oscillations(isnan(oscillations.dropID),:)=[];
        
        oscillations_peaks=oscillations.peakTime;
        oscillations_dropids=oscillations.dropID;
        oscillations_pperiods=oscillations.periodPeak;
        oscillations_cycleids=oscillations.cycleID;
        
        for o=1:size(oscillations,1)
            % info of an oscillation of interest
            tmp_dropID=oscillations_dropids(o);
            tmp_peak=oscillations_peaks(o);
            tmp_pperiods=oscillations_pperiods(o);
            tmp_cycleid=oscillations_cycleids(o);
            
            % info of tracking result of droplet
            track=tracks(tracks_ids==tmp_dropID).Feat;
            track_start=tracks(tracks_ids==tmp_dropID).Start;
            track_Nuc=tracks(tracks_ids==tmp_dropID).Nuc;
            track_Len=tracks(tracks_ids==tmp_dropID).Len;
            track_nucarea=track.Nuc_area;
            track_cytoarea=track.Cyto_area;
            track_time=track.t;
            
            % adding data.Nuc for droplet level nucleus existence
            if track_Nuc > track_Len*thresh/100
                oscillations.Nuc(o)=1;
            else
                oscillations.Nuc(o)=0;
            end
            
            % adding data.duration, data.seg_frame, data.seg_ratio for 
            % single oscillation level nuclear detection
            if tmp_cycleid==1 && isnan(tmp_pperiods)
                oscillations.duration(o)=tmp_peak;
                tmp_p2_track_index=find(track_time==tmp_peak);
                tmp_track_nucarea=track_nucarea(1:tmp_p2_track_index);
                tmp_track_cytoarea=track_cytoarea(1:tmp_p2_track_index);
                
            else
                oscillations.duration(o)=tmp_pperiods;
                tmp_p1_track_index=find(track_time==tmp_peak-tmp_pperiods);
                tmp_p2_track_index=find(track_time==tmp_peak);
                tmp_track_nucarea=track_nucarea(tmp_p1_track_index:tmp_p2_track_index);
                tmp_track_cytoarea=track_cytoarea(tmp_p1_track_index:tmp_p2_track_index);
            end
            cmp_res=tmp_track_nucarea < tmp_track_cytoarea; % compare
            oscillations.seg_frame(o)=sum(cmp_res);            
            oscillations.seg_ratio(o)=oscillations.seg_frame(o)/oscillations.duration(o);
        end
        
    end
   
    %update data
    oscillation_table(positionIndex).data=oscillations;
    positionIndex=positionIndex+1;
end

end

