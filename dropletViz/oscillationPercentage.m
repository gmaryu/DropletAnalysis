function [osci_percentage_array, osci_percentage_nuc, osci_percentage_cyto]=oscillationPercentage(position, sumTrackTable, oscillation_table, option,group)
%{
    Calculation of oscillation percentage
    position: a list of position of you need to plot
    sumTrackTable: tracking result of all detected doroplets
    oscillation_table: oscillation info summary of all droplets
    option: plotting option, assign integer below 
        1 : for all droplets (one barplot)
        2 : for all, nucleus, cytoplasm (three barplots)
    group: a list for color assignment to plot different conditions
    
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
osci_percentage_array=[]; % for all droplets
osci_percentage_nuc=[];   % nucleus droplets
osci_percentage_cyto=[];  % cytoplasm droplets
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
        osci_percentage_array=[osci_percentage_array;tmp_percentage];
        
        if option == 2
            tmp_nuc_detected=extractfield(tmp_tracks, 'Nuc');
            tmp_tracked_len=extractfield(tmp_tracks, 'Len');
            thres=tmp_tracked_len.*0.3; % more than 30%
            
            tmp_tracks_nuc=tmp_tracks(tmp_nuc_detected>thres);
            tmp_tracks_cyto=tmp_tracks(tmp_nuc_detected<=thres);
            
            if ~isempty(tmp_tracks_nuc)
                tmp_tracks_nuc_dropletid=extractfield(tmp_tracks_nuc, 'id');
                len_osci_nuc=length(intersect(tmp_tracks_nuc_dropletid, unique(tmp_oscillationTbl.data.dropID)));
                per_osci_nuc=len_osci_nuc/length(tmp_tracks_nuc_dropletid);
                osci_percentage_nuc=[osci_percentage_nuc; per_osci_nuc];
            else
                osci_percentage_nuc=[osci_percentage_nuc; NaN];
            end
            
            tmp_tracks_cyto_dropletid=extractfield(tmp_tracks_cyto, 'id');
            len_osci_cyto=length(intersect(tmp_tracks_cyto_dropletid, unique(tmp_oscillationTbl.data.dropID)));
            per_osci_cyto=len_osci_cyto/length(tmp_tracks_cyto_dropletid);
            osci_percentage_cyto=[osci_percentage_cyto; per_osci_cyto];
        end
        tlabel=horzcat(tlabel, int2str(allpositions(i)));
    end
end



figure();
b1=bar(osci_percentage_array);
for c=1:length(group)
    cc=color_code{group(c)};
    b1.FaceColor = 'flat';
    b1.CData(c,:) = cc;
end
hold on
set(gca,'xticklabel',tlabel);
title('All droplets')
ylabel('Percentage')
xlabel('Position')
hold off

if option == 2
    figure();
    b2=bar(osci_percentage_nuc);
    for c=1:length(group)
        cc=color_code{group(c)};
        b2.FaceColor = 'flat';
        b2.CData(c,:) = cc;
    end
    hold on
    set(gca,'xticklabel',tlabel);
    title('Nuclear droplets');
    ylabel('Percentage');
    xlabel('Position');
    hold off
    
    figure();
    b3=bar(osci_percentage_cyto);
    for c=1:length(group)
        cc=color_code{group(c)};
        b3.FaceColor = 'flat';
        b3.CData(c,:) = cc;
    end
    hold on
    set(gca,'xticklabel',tlabel)
    title('Cytoplasm droplets');
    ylabel('Percentage');
    xlabel('Position');
    hold off

end