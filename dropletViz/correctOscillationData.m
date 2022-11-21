function correctOscillationData(dataPath, savePath, position, targetPos, sumTrackTable)
%{
Description
    manual peak detection in a designated position

Args
    dataPath: original directory of peak detection result
    savePath: destination directory for saving this result
    position: a list of all positions performed droplet segmentation
    targetPos: a position for this analysis (single position accepted)
    sumTrackTable: cleaned tracking result
%}

%% peakCollection

% Load data
%{
dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20210203_Aphidicolin_rb\peakSelection';
savePath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20210203_Aphidicolin_rb'; %save fig path
pos=[6];
%}
osci_table=readOscillationData(dataPath, targetPos);
osci_table2=osci_table.data; % return value, to be updated table
tracks=sumTrackTable{position==osci_table.position};
droplets=extractfield(tracks,'id');

%% Main
%del_tbl=zeros(size(osci_table2,1),1);
targetDID=unique(osci_table2.dropID); % target droplet ID
%targetDID=[7,40,62]; % target droplet ID
%targetDID=[]; % target droplet ID
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
    tblIDX=find(osci_table.data.dropID==targetDID(i)); % get index of target dropletID in origianl table for data delition from updated table
    
    %disp(tblIDX);
    peakTime=osci_table.data.peakTime(tblIDX)+min(t)-1; % get peak timepoints (adjust to image acquisition frame number)
    peakValue=osci_table.data.peakValue(tblIDX); % get peak values
    tmpdTable=osci_table.data(tblIDX,:); % extract target droplet oscillation data 
    
    %% figure plot
    figure(1);
    hold on
    plot(t, sigValue); % time lapse imaging data from "tracks"
    scatter(peakTime, peakValue); % previously detected peak information
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
            [vv,idx]=max(sigValue(neighbor_t-min(t)+1));
            tt=neighbor_t(idx);
            if ~isempty(find(peakTime==tt, 1))
                warning('peak already in list')
                %return;
            end
            
            % for plot
            peakTime=[peakTime;tt];
            peakValue=[peakValue;vv];
            figure(1)
            clf('reset')
            hold on
            plot(t, sigValue);
            scatter(tt,vv,60,'*');
            scatter(peakTime, peakValue);
            hold off
            
            % for table update
            del_tbl(tblIDX)=1; % delete all row from orignal position
            tmp=NaN*ones(1,size(tmpdTable,2)); % new row for selected peak
            tmp(1)=tt; %peakTime
            tmp(2)=vv; %peakValue
            tmpdTable=[tmpdTable; array2table(tmp,'VariableNames',osci_table2.Properties.VariableNames)]; % add row to target droplet
            
        elseif button==3 %Delete peak: right click
            e=1;
            % search delete peak from selected point
            [~,idx]=min(abs(peakTime-x)); % get index of closest peak
            %tt=peakTime(idx);
            %vv=peakValue(idx);
            
            % for plot
            figure(1)
            clf('reset')
            peakTime(idx)=[]; % delete selected peak data from list 
            peakValue(idx)=[]; % delete selected peak data from list
            hold on
            plot(t, sigValue);
            scatter(peakTime,peakValue,60)
            hold off
            
            % for table update
            del_tbl(tblIDX)=1; % delete all row from orignal position
            %del_tbl(tblIDX(idx))=1;
            tmpdTable(idx,:)=[];
            
        elseif button=='d' % remove all peaks
            e=1;
            % for plot
            peakTime(:)=[];
            peakValue(:)=[];
            figure(1)
            clf('reset')
            hold on
            plot(t, sigValue);
            scatter(peakTime,peakValue,60)
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
    
    if e==0
        tmpdTable(:,:)=[]
    end
    % modify table
    if e==1 && ~isempty(tmpdTable)
        tmpdTable=sortrows(tmpdTable,'peakTime');
        tmpdTable.cycleID=[1:size(tmpdTable,1)]';
        tmpdTable.dropID=mode(tmpdTable.dropID)*ones(size(tmpdTable,1),1);
        tmpdTable.posIdx=mode(tmpdTable.posIdx)*ones(size(tmpdTable,1),1);
        tmpdTable.position=mode(tmpdTable.position)*ones(size(tmpdTable,1),1);
        tmpdTable.frameStep=mode(tmpdTable.frameStep)*ones(size(tmpdTable,1),1);
        tmpdTable.periodPeak=[NaN;tmpdTable.peakTime(2:end)-tmpdTable.peakTime(1:end-1)];
        add_table=[add_table;tmpdTable];
    end
    
end
osci_table2(logical(del_tbl),:)=[];
osci_table2=[osci_table2;add_table];


%% save mat file
formatOut=30;
imPath=fullfile(savePath, 'PeakSelectionCorrected', ['Pos',num2str(targetPos)]);
mkdir(imPath);
save(fullfile(imPath,['stats_summary_',datestr(now,formatOut),'.mat']),'tracks','osci_table2','-v7.3');
end