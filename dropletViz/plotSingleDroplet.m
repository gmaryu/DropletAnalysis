function plotSingleDroplet(sumTrackTable, allpositions, position, cell_id, label, ax, interval)
%{
 Description: plot a single plot fot designated droplet
 Params:
    sumTrackTable: cleaned tracking result
    allposition: a list of all positions to make a sumTrackTable
    position: a paticular position number for a droplet
    cell_id: a cell_id 
    label: a cell structure of column names
        e.g.: {FC_ration_whole, CFPstdInt}
    ax: a list of axis "1" is for left, "2" is for right
        e.g.: [1,2]
    interval: time interval (min) between time posinits
%}
%% 
if length(sumTrackTable) ~= length(allpositions)
    disp("length of sumTrackTable and all positions is different");
    return;
end
selected_position=ismember(allpositions, position);

%% data extraction
for s=1:length(selected_position)
   if selected_position(s)==1
       tmpTracks=sumTrackTable{s};
       tmp_ids=extractfield(tmpTracks, 'id');
       [v, selected_droplet_index]=max(ismember(tmp_ids,cell_id));
       if v == 1
           tmp_data=tmpTracks(selected_droplet_index).Feat;
       else
           disp('this cell_id is not found.')
           break;
       end
       
   end
end

%% plot
if length(label)~=length(ax)
    disp('length of label and ax is mismatch.')
    return
end
figure();hold on;
for l=1:length(label)
    xdata=tmp_data.t*interval;
    ydata=eval(['tmp_data.',label{l},';']);
    if ax(l)==1
        yyaxis left;
    else
        yyaxis right;
    end
    plot(xdata,ydata);
end
hold off;


