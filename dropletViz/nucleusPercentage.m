function nucleusPercentage(positions, qpositions, sumTrackTable, nucdetected)
%{
    Calculation of a percentage of nucleus droplets
    positions: a list of position of you need to plot
    qpositions: a list of positions where you quantified.like "genPos"
    sumTrackTable: tracking result of all detected doroplets
    nucdetected: threshold for nucleus droplet definition. Percentage of
    nuclues detected frame number.
    
    e.g.
        nucleusPercentage([10,11,12,13,14,15], genPos, sumTrackTable)
%}

%%
t=nucdetected/100;
%%
percentage_nuc=[];   % nucleus droplets
percentage_cyto=[];  % cytoplasm droplets

% selected position will be calculate
selected_position=ismember(qpositions, positions);

for i=1:length(selected_position)
    if selected_position(i)==1
        % all droplets
        tmp_tracks=sumTrackTable{i};
        num_tracked_droplets=length(tmp_tracks);
        tmp_nuc_detected=extractfield(tmp_tracks, 'Nuc');
        tmp_tracked_len=extractfield(tmp_tracks, 'Len');
        thres=tmp_tracked_len.*t; % more than 30%
        
        tmp_tracks_nuc=tmp_tracks(tmp_nuc_detected>thres);
        tmp_tracks_cyto=tmp_tracks(tmp_nuc_detected<=thres);
        
        if ~isempty(tmp_tracks_nuc)
            per_tracked_nuc=length(tmp_tracks_nuc)/num_tracked_droplets;
            percentage_nuc=[percentage_nuc; per_tracked_nuc];
            per_tracked_cyto=length(tmp_tracks_cyto)/num_tracked_droplets;
            percentage_cyto=[percentage_cyto; per_tracked_cyto];
        else
            disp(['No nucleus droplets in Position', num2str(qpositions(i))]);
            percentage_nuc=[percentage_nuc];
            percentage_cyto=[percentage_cyto];
        end
    end
end

figure();
pie([mean(percentage_nuc), mean(percentage_cyto)]);
legend({'Nucleus','Cytoplasm'})
title('Percentage of nucleus droplets')
