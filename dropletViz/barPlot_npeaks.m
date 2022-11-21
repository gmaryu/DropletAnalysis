function [npeaks_data, npeaks_stdev, b]=barPlot_npeaks(position, oscillation_table, group)
%{
Description
    plot a bar plot for the number of peaks (including 1st peak) in each
    position
Args
    pos: position what you plot
    osci_table: oscillation data table
    group: group of experiment condition. assign as integer
        e.g. [1,1,1,2,2,2,3,3,3,4,4,4,4,4]
%}
%%
color_code={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
            [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
npeaks_data=[];
npeaks_stdev=[];
tlabel={};
allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, position);
for s=1:length(selected_position)
    if selected_position(s)==1
        PosN_npeaks=[];
        PosN_osci=oscillation_table(s).data;
        %eval(['PosN_osci=Pos',int2str(pos(posIdx)),'_osci;']);
        PosN_unq_dropletID=unique(PosN_osci.dropID);
        for i=1:length(PosN_unq_dropletID)
            npeaks=find(PosN_osci.dropID==PosN_unq_dropletID(i)); % number of peaks
            PosN_npeaks=[PosN_npeaks length(npeaks)];
        end
        eval(['Pos',int2str(allpositions(s)),'_npeaks=PosN_npeaks;']);
        npeaks_data=[npeaks_data, mean(PosN_npeaks)];
        npeaks_stdev=[npeaks_stdev, std(PosN_npeaks)];
        tlabel=horzcat(tlabel, int2str(allpositions(s)));
    end
end
figure();
b=bar(1:length(npeaks_data), npeaks_data);
%%
for c=1:length(group)
    cc=color_code{group(c)};
    b.FaceColor = 'flat';
    b.CData(c,:) = cc;
end
hold on
er = errorbar(1:length(npeaks_data),npeaks_data,npeaks_stdev);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'xticklabel',tlabel)
ylabel('Cycle Number')
xlabel('Position')
%hold off
%savefig('barPlot_npeaks.fig')
end