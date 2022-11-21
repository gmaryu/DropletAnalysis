function [osci_npeaks_array, osci_npeaks_stdev]=plotCycleNumber(position, sumTrackTable, oscillation_table, option, group)
%{
    plot peak number in each position
    position: a list of position of you need to plot
    sumTrackTable: tracking result of all detected doroplets
    oscillation_table: oscillation info summary of all droplets
    option: plotting option, assign integer below 
        1 : for all droplets (one barplot)
        2 : for all, nucleus, cytoplasm (three barplots)
    group: a list for color assignment to plot different conditions. 
           This option doesn't work on option '2'.
    
    e.g.
        [osci_percentage_array]=oscillationPercentage(genPos, sumTrackTable, oscillation_table, 2, [1,1,2,2,3,3,3,3,3,3]);
%}
%% check if the selected option is applicable
if ~isfield(sumTrackTable{1}, 'Nuc') && option==2
    disp('No nucleus segmentation results');
    return
end
%%
color_code={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
            [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
%%
t =0.2;
osci_npeaks_array=[]; % for all droplets
osci_npeaks_stdev=[];
osci_npeaks_nuc=[];   % nucleus droplets
osci_npeaks_nuc_stdev=[];
osci_npeaks_cyto=[];  % cytoplasm droplets
osci_npeaks_cyto_stdev=[];
tlabel={};

allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, position);

if length(sumTrackTable)~=length(oscillation_table)
    disp('Number of position of Tracking results and  is mismatch.');
end

for i=1:length(sumTrackTable)
    if selected_position(i)==1
        % all droplets
        tmp_tracks=sumTrackTable{i};
        tmp_oscillationTbl=oscillation_table(i);
        num_tracked_droplets=length(tmp_tracks);
        num_osci_droplets=length(unique(tmp_oscillationTbl.data.dropID));
        tmp_percentage=num_osci_droplets/num_tracked_droplets;
        
        % all droplets
        PosN_npeaks=[];
        PosN_osci=tmp_oscillationTbl.data;
        %eval(['PosN_osci=Pos',int2str(pos(posIdx)),'_osci;']);
        PosN_unq_dropletID=unique(PosN_osci.dropID);
        for j=1:length(PosN_unq_dropletID)
            npeaks=find(PosN_osci.dropID==PosN_unq_dropletID(j)); % number of peaks
            PosN_npeaks=[PosN_npeaks length(npeaks)];
        end
        osci_npeaks_array=[osci_npeaks_array, mean(PosN_npeaks)];
        osci_npeaks_stdev=[osci_npeaks_stdev, std(PosN_npeaks)];
            
        if option == 2
            tmp_nuc_detected=extractfield(tmp_tracks, 'Nuc');
            tmp_tracked_len=extractfield(tmp_tracks, 'Len');
            thres=tmp_tracked_len.*t; % more than 30%
            
            tmp_tracks_nuc=tmp_tracks(tmp_nuc_detected>thres);
            tmp_tracks_cyto=tmp_tracks(tmp_nuc_detected<=thres);
            
            if ~isempty(tmp_tracks_nuc)
                tmp_tracks_nuc_dropletid=extractfield(tmp_tracks_nuc, 'id');
                for k =1:length(tmp_tracks_nuc_dropletid)
                    npeaks=find(PosN_osci.dropID==tmp_tracks_nuc_dropletid(k)); % number of peaks
                    PosN_npeaks=[PosN_npeaks length(npeaks)];
                end
                osci_npeaks_nuc=[osci_npeaks_nuc, mean(PosN_npeaks)];
                osci_npeaks_nuc_stdev=[osci_npeaks_nuc_stdev, std(PosN_npeaks)];
                
            else
                disp(['No nucleus droplets in Position', num2str(allpositions(i))]);
                osci_npeaks_nuc=[osci_npeaks_nuc, NaN];
                osci_npeaks_nuc_stdev=[osci_npeaks_nuc_stdev, NaN];
            end
            
            tmp_tracks_cyto_dropletid=extractfield(tmp_tracks_cyto, 'id');
            for k =1:length(tmp_tracks_cyto_dropletid)
                npeaks=find(PosN_osci.dropID==tmp_tracks_cyto_dropletid(k)); % number of peaks
                PosN_npeaks=[PosN_npeaks length(npeaks)];
            end
            osci_npeaks_cyto=[osci_npeaks_cyto, mean(PosN_npeaks)];
            osci_npeaks_cyto_stdev=[osci_npeaks_cyto_stdev, std(PosN_npeaks)];
        end
        tlabel=horzcat(tlabel, int2str(allpositions(i)));
    end
end

if option == 1
    figure();
    b=bar(1:length(osci_npeaks_array), osci_npeaks_array);
    for c=1:length(group)
        cc=color_code{group(c)};
        b.FaceColor = 'flat';
        b.CData(c,:) = cc;
    end
    hold on
    er = errorbar(1:length(osci_npeaks_array),osci_npeaks_array,osci_npeaks_stdev);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    set(gca,'xticklabel',tlabel)
    ylabel('Cycle Number')
    xlabel('Position')
    %hold off
elseif option == 2
    figure();
    nepaks_array=cat(1, osci_npeaks_nuc, osci_npeaks_cyto, osci_npeaks_array);
    stdev_array=cat(1, osci_npeaks_nuc_stdev, osci_npeaks_cyto_stdev, osci_npeaks_stdev);
    b=bar(1:length(osci_npeaks_array), nepaks_array);
    set(gca,'xticklabel',tlabel)
    legend('Nucleus','Cytoplasm','All')
    ylabel('Cycle Number')
    xlabel('Position')    
end

end