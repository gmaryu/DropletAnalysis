
function [binarysignalmat, sortdata, sorteddata]=rasterBinaryMat(oscillation_table, position, timepoints, droplets_cell, sortcolumn)
%{
oscillation_table=oscillation_table;
position=[1,2,3,4];
timepoints=375;
droplets_cell: cell of "list of droplets id" or 0. 0 is for all droplets. 
sortcolumn='DyeIntensity' or ''

ex.
rasterBinaryMat(osci_table, pos, 400, {'all','all','all','all','all'}, '');
%}
%%
%posIdx=1;
%npeaks_data=[]; %
%npeaks_stdev=[];
binarysignalmat=[]; % 1: peak, 0: other
sortdata=[]; % sorted dye intensity vector

if length(position) ~= length(droplets_cell)
    disp("Length of position and droplets_cell don't match.");
    return;
else
    figure();
    hold on
    fidx=[]; % field index
    pnum=[]; % position number
    for n=1:length(oscillation_table)
        fidx=[fidx n]; % field index
        pnum=[pnum oscillation_table(n).position];
    end
    
    for i=1:length(position)
        f=pnum==position(i);% field index of target position
        PosN_npeaks=[];
        PosN_osci=oscillation_table(fidx(f)).data; % load oscillations data
        if droplets_cell{i} ~= 0
            PosN_unq_dropletID = cell2mat(droplets_cell(i));
        else
            PosN_unq_dropletID=unique(PosN_osci.dropID); % get all doplets id that have oscillations
        end
        for i=1:length(PosN_unq_dropletID)
            dropletID=PosN_unq_dropletID(i); %target id
            tmp=zeros(1,timepoints); % make blank vector for peak information
            npeaks=PosN_osci.dropID==dropletID; % logical vector that indicates peak position in table
            tpdata=PosN_osci.peakTime(npeaks); % vector of timepoints
            if ~isempty(sortcolumn)
                eval(['cdata=PosN_osci.',sortcolumn,';']); % column data
                sortdata=[sortdata, min(cdata(npeaks))];
            end
            tmp(tpdata)=1;
            binarysignalmat=vertcat(binarysignalmat,tmp);
        end
    end
    [sorteddata,sortedorder]=sort(sortdata);
    if ~isempty(sortcolumn)
        binarysignalmat=binarysignalmat(sortedorder,:);
    end
    binarysignalmat=logical(binarysignalmat);
    LineFormatHorz.LineWidth = 3;
    LineFormatHorz.Color = 'b';
    LineFormatVert.LineWidth = 5;
    MarkerFormat.MarkerSize = 12;
    MarkerFormat.Marker = '*';
    f=plotSpikeRaster(binarysignalmat,'PlotType','vertline','LineFormat',LineFormatVert,'VertSpikePosition', -0.5);
    xlabel('Time frame')
    ylabel('Indivisual droplets')
end
end