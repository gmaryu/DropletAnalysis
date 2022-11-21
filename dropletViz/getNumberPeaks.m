function npeaks_list=getNumberPeaks(npeaks_list, tmp_oscillation_table)


tmp_unq_dropletID=unique(tmp_oscillation_table.data.dropID);
for i=1:length(tmp_unq_dropletID)
    npeaks=find(tmp_oscillation_table.data.dropID==tmp_unq_dropletID(i)); % number of peaks
    npeaks_list=[npeaks_list length(npeaks)];
end

