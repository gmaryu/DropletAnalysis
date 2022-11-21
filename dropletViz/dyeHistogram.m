function h = dyeHistogram(pos_groups, sumTrackTable, allPosition, legend_groups)
%{
Description
    plot a histogram of dye intensity

Args
    pos_groups: a cell list of positions where you want to plot
        e.g. pos_groups={[6,7],[8,9],[10,11,12,13,14,15]};
    sumTrackTable: tracking results that include dye intensity data
    allposition: list of position. This must be the same size as sumTrackTable
    legend_groups: a cell list of legends
        e.g. legend_groups={'A','B','C'};
%}

upper_list=[];
lower_list=[];
%%
if length(allPosition) ~= length(sumTrackTable)
    disp('Params have different length');
end

%%
for g = 1:length(pos_groups)
    tmp_group=pos_groups{g};
    dyeInt=[];
    for p = 1:length(tmp_group)
        selected_position=ismember(allPosition, tmp_group(p));
        for tt=1:length(sumTrackTable)
            if selected_position(tt) == 1
                tracks=sumTrackTable{tt};
                for t=1:length(tracks)
                    dye=tracks(t).Feat.BFP_whole(1);
                    dyeInt=[dyeInt, dye];
                end
            end
        end
        
    end
    upper_list=[upper_list, prctile(dyeInt,98)];
    lower_list=[lower_list, prctile(dyeInt,2)];
    h=histogram(dyeInt);
    h.BinWidth=200;
    hold on
end

xlim([min(lower_list),max(upper_list)]);
legend(legend_groups);
title('Distrubution of dye intensity')
ylabel('Droplet count')
xlabel('Dye Intensity')
hold off
end