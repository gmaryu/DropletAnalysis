function correctTroughData(peakcorrected_data, savePath, targetPos)
%{
Description
    manual trough detection in a designated position

Args
    peakcorrected_data: 
    savePath: destination directory for saving this result
    targetPos: a position for this analysis (single position accepted)
%}

%% troughCollection

% Load data
%{
dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20210203_Aphidicolin_rb\peakSelection';
savePath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20210203_Aphidicolin_rb'; %save fig path
pos=[6];
%}
%osci_table=readOscillationData(dataPath, targetPos);
%osci_table2=osci_table.data; % return value, to be updated table
osci_table2=peakcorrected_data.osci_table2;
%tracks=sumTrackTable{position==osci_table.position};
tracks=peakcorrected_data.tracks;
droplets=extractfield(tracks,'id');

%% Main
%del_tbl=zeros(size(osci_table2,1),1);
targetDID=unique(osci_table2.dropID); % target droplet ID
targetDID=targetDID(~isnan(targetDID));
%targetDID=[7,40,62]; % target droplet ID

add_table=table(); % all edited droplet informatoin 
del_tbl=zeros(size(osci_table2,1),1); % use as logical; for delete all edited droplet information

for i=1:length(targetDID) % for loop by peak detected droplet in automation
    %%
    disp(targetDID(i)) % previously peak detected dropeltID
    dropletIDX=droplets==targetDID(i); % get index at "tracks"
    %disp(dropletIDX);
    track=tracks(dropletIDX).Feat; % extract quantified data of target dropelt
    sigValue=track.FC_ratio_whole; %should be modified later
    t=track.t; % timepoint data (frame)
    tblIDX=find(osci_table2.dropID==targetDID(i)); % get index of target dropletID in origianl table for data delition from updated table
    
    %disp(tblIDX);
    peakTime=osci_table2.peakTime(tblIDX)+min(t)-1; % get peak timepoints (adjust to image acquisition frame number)
    peakValue=osci_table2.peakValue(tblIDX); % get peak values
    troughTime=osci_table2.troughTime(tblIDX)+min(t)-1; % get peak timepoints (adjust to image acquisition frame number)
    troughValue=osci_table2.troughValue(tblIDX); % get peak values
    tmpdTable=osci_table2(tblIDX,:); % extract target droplet oscillation data 
    
    if length(peakTime)>1
        %% figure plot
        figure(1);
        hold on
        plot(t, sigValue); % time lapse imaging data from "tracks"
        scatter(peakTime, peakValue,'v');
        scatter(troughTime, troughValue); % previously detected peak information
        hold off
        %% data correction
        [x,y,button] = ginput(1);
        e=0;
        while button~='f'
            if button==1 %Add peak: left click
                % search max value point from selected point
                e=1;
                r1=floor(x-5);
                r2=floor(x+5);
                if r1 <= 0
                    neighbor_t=(1:floor(x+5)); % real-timescale
                elseif r2> max(t)
                    neighbor_t=(floor(x-5):max(t)); % real-timescale
                else
                    neighbor_t=(floor(x-5):floor(x+5)); % real-timescale
                    %disp(min(t)-1) % 1st frame -1
                end
                [vv,idx]=min(sigValue(neighbor_t-min(t)+1));
                tt=neighbor_t(idx);
                if ~isempty(find(troughTime==tt, 1))
                    warning('trough already in list')
                    %return;
                end
                
                % for plot
                troughTime=[troughTime;tt];
                troughValue=[troughValue;vv];
                figure(1)
                clf('reset')
                hold on
                plot(t, sigValue);
                scatter(tt,vv,60,'*');
                scatter(troughTime, troughValue);
                hold off
                
                % for table update
                del_tbl(tblIDX)=1; % delete all row from orignal position
                pp=tmpdTable.peakTime>tt;%logical
                tind=find(pp==1, 1);
                tmpdTable(tind,'troughTime')={tt};
                tmpdTable(tind,'troughValue')={vv};
                %disp(tmpdTable);
                %tmpdTable=[tmpdTable; array2table(tmp,'VariableNames',osci_table2.Properties.VariableNames)]; % add row to target droplet
                
            elseif button==3 %Delete peak: right click
                e=1;
                % search delete peak from selected point
                [~,idx]=min(abs(troughTime-x)); % get index of closest peak
                %tt=peakTime(idx);
                %vv=peakValue(idx);
                
                % for plot
                figure(1)
                clf('reset')
                troughTime(idx)=NaN; % delete selected peak data from list
                troughValue(idx)=NaN; % delete selected peak data from list
                hold on
                plot(t, sigValue);
                scatter(troughTime,troughValue,60)
                hold off
                
                % for table update
                del_tbl(tblIDX)=1; % delete all row from orignal position
                %del_tbl(tblIDX(idx))=1;
                tmpdTable(idx,'troughTime')={NaN};
                tmpdTable(idx,'troughValue')={NaN};
                tmpdTable(idx,'ampLeft')={NaN};
                tmpdTable(idx,'ampRight')={NaN};
                tmpdTable(idx,'periodTrough')={NaN};
                disp(tmpdTable);
                
            elseif button=='d' % remove all peaks
                e=1;
                % for plot
                troughTime(:)=[];
                troughValue(:)=[];
                figure(1)
                clf('reset')
                hold on
                plot(t, sigValue);
                scatter(troughTime,troughValue,60)
                hold off
                
                % for table update
                del_tbl(tblIDX)=1;
                tmpdTable(:,:)=[];
            else
                disp('this action is invalid.')
                
            end
            
            [x,y,button] = ginput(1);
        end
        close gcf;
    else
        continue;
    end
    if e==0
        tmpdTable(:,:)=[];
    end
    % modify table
    if e==1 && ~isempty(tmpdTable)
        tmpdTable.ampLeft=[NaN;tmpdTable.peakValue(2:end)-tmpdTable.troughValue(2:end)];
        tmpdTable.ampRight=[NaN;tmpdTable.peakValue(1:end-1)-tmpdTable.troughValue(2:end)];
        tmpdTable.periodTrough=[tmpdTable.troughTime(2:end)-tmpdTable.troughTime(1:end-1);NaN];
        add_table=[add_table;tmpdTable]
        disp(tmpdTable);
    end
    
end
osci_table2(logical(del_tbl),:)=[];
osci_table2=[osci_table2;add_table];


%% save mat file
formatOut=30;
imPath=fullfile(savePath, 'TroughSelectionCorrected', ['Pos',num2str(targetPos)]);
mkdir(imPath);
save(fullfile(imPath,['stats_summary_',datestr(now,formatOut),'.mat']),'tracks','osci_table2','-v7.3');
end