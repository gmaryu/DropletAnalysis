% use this script after data cleaning
%% params
% save dir of analysis results
dataPath='PATH_TO_RESULT_DIRECTORY'; % store directory

% name of sub-directory
plotType={'peakSelection'};

% position for analysis
position=genPos;   % all position to be processedposition=genPos;   % all position to be processed
posIdx=1;

% channel for peak detection
colNames={'FC_ratio_whole'};
refNames={'FRET_whole','CFP_whole'}; % all position to be processed
sigScale=1;

mannualEdit=0;          % whether do mannual edition
oldF=0;                 % old naming format
loadSaveFile=0;         % load saved file?(same as dataPath)

% data save option
isSaveTrackFig=0;       %save editing data to fig?
isSaveTrackPng=1;       %save editing data to png?
isPlotStats=1;          %plot statistics figure
formatOut=30;

% params for peak detection
frameStep=3;        % min,frame rate, aquisition interval
ampThresh=[-2,5];   % binding range from smooth data to real data
adaptWindSize=3;    % binding window size from smoothed data to real data
timeThresh=1;      % Minimum time point for a track to be eligible for analysis
senThresh=0.01;     % min prom of envolope difference 

%prom: prominence?
%promSens=5;       % added to envolope difference
promSens=0.06; 
smoothThresh=20;    % window size for outlier detection
peakMinDist=5;     % minimum peak distance
peakMaxDist=400;     % /2 for generation of envolope

isClear=0;          % flag for clear and position
filtThresh=[2,0.5,5]; % max linear fitting(time vs amp) resid; max left amp/right amp; max prom/(max previous prom);

%robust linear fitting residual; left right amplitude ratio; max before/min
%amplitude
lbank=0.5/peakMinDist;
lpFilt = designfilt('lowpassiir','PassbandFrequency',lbank, ...
    'StopbandFrequency',lbank+0.05,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65);% low pass filter reserved
% baseFilt = designfilt('lowpassiir','PassbandFrequency',0.005, ...
%     'StopbandFrequency',0.025,'PassbandRipple',1, ...
%     'StopbandAttenuation',10);
tempSig=signalClassPT(senThresh,promSens,smoothThresh,peakMinDist,filtThresh,...
    lpFilt,peakMaxDist,adaptWindSize,ampThresh);   
trackCell=[];
drop=[];
osci=[];
raw=[];
peakThresh=1;
skipFlag=0;
osciDat=[];
dropDat=[];
rawDat=[];

%% load saved file
if loadSaveFile
    saveFileName=ls([statsPath,'\stats_track_*.mat']);
    %saveFileName=ls([dataPath,'Stat_Pos',num2str(position(1)),'_*.mat']);
    if isempty (saveFileName)
        warning('no save file found')
        return;
    else
        saveFileName=saveFileName(end,:);
        disp(['loading ',saveFileName])
        load(fullfile(statsPath,saveFileName),'osci','drop','raw')
        saveFile=load(fullfile(statsPath,saveFileName),'tempSig','ii','posIdx','isClear');
        clearFlag=saveFile.isClear;
    end
end

%% main loop by the position
while posIdx <= length(position)
    disp('aaaa');
    % overwrite position index if loadSaveFile was true
    if posIdx==1 && loadSaveFile==1
        posIdx=saveFile.posIdx;
    end
    
    % load data
    for n=1:size(posList{1},2)
        if strcmp(['Pos',num2str(position(posIdx))],posList{1}{n})
            relPos=n;
            disp(relPos);
            % for save results
            imPath=fullfile(dataPath, plotType{1}, posList{1}{n});
            disp(imPath);
            mkdir(imPath)
            rawPath=fullfile(imPath,'raw');
            mkdir(rawPath);
            statsPath=fullfile(imPath,'stats');
            disp(statsPath);
            mkdir(statsPath);
            break;
        end
    end
    
    % load tracking result
    tracks=sumTrackTable{n};
 
    % loop by droplets
    ii=1;
    while ii <= length(tracks)
        disp(['ii: ',num2str(ii)]); % droplet id
        
        if ii==1&&loadSaveFile==1&&~clearFlag
            % load saved matfile 
            disp(saveFile.ii); % last droplet number in saved mat file
            ii=saveFile.ii; % overwrite 
            tempSig=saveFile.tempSig;
            %disp(['1,load,clear ',num2str(ii)]);
        else
            if ii==1 && loadSaveFile==1
                ii=saveFile.ii;
                %disp(ii);
                %disp(['1,load  ',num2str(ii)]);
            end
            
            eval(['aveSig=tracks(',num2str(ii),').Feat.',colNames{1},';']);
            eval(['madSig=tracks(',num2str(ii),').Feat.',colNames{1},';']);
            aveSig=aveSig*sigScale; %
            madSig=madSig*sigScale; % unused? 
            
            areaSig=tracks(ii).Feat.area;
            xcoorSig=tracks(ii).Feat.xcoord;
            ycoorSig=tracks(ii).Feat.ycoord;
            tSig=tracks(ii).Feat.t;
            
            tempSig.setRef(eval(['tracks(',num2str(ii),').Feat.',refNames{1},';']),...
                           eval(['tracks(',num2str(ii),').Feat.',refNames{2},';'])); % unused?
            tempSig.setOri(aveSig); % target signal. original and processed signal are preped. find peak trough with raw data.
            tempSig.feature=table(tSig,madSig,areaSig,xcoorSig,ycoorSig);
            tempSig.centerSig();
            tempSig.fillOut(); % filloutliers with smooth thresh; doesn't have much effect
            tempSig.genEnv(); 
            tempSig.genPeak();
            tempSig.attachPeak();
            %disp(length(tempSig.peakClass.time));
            
            if length(tempSig.peakClass.time)<peakThresh
                ii=ii+1;
                tempSig.clearObj();
                disp('below peakThresh');
                continue;
            end
        end
        tempSig.showPeak();
        %% selection peaks
        if mannualEdit==1
            [x,y,button] = ginput(1);
            while button~='2'
                if button==1 %left click
                    tempSig.userAddPeak(x);
                elseif button==3 %right click
                    tempSig.userDelPeak(x);
                elseif button==2 %
                    tempSig.peakClass=tempSig.backClass;
                    tempSig.showPeak();
                elseif button=='s'
                    save([statsPath,'/stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                    disp('saving finished');
                elseif button=='r'
                    % does this work?
                    tempSig.lot
                    Peak();
                elseif button=='q' %exit
                    return
                elseif button=='e' %skip
                    tempSig.clearObj();
                    skipFlag=1;
                    break;
                elseif button=='c'
                    % does this work?
                    isClear=1;
                    save([statsPath,'/stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                    disp('saving finished for previous droplets');
                end
                [x,y,button] = ginput(1);
            end
            if skipFlag
                skipFlag=0;
                ii=ii+1;
                disp(['skip_peak: ', num2str(ii)]);
                continue;
            end
        end
        
        
        %% trough selection
        tempSig.genTrough();
        tempSig.showTrough();
        
        if mannualEdit==1
            [x,y,button] = ginput(1);
            while button~='2'
                if button==1
                    tempSig.userAddTrough(x);
                elseif button==3
                    tempSig.userDelTrough(x);
                elseif button==2
                    tempSig.troughClass=tempSig.backClass;
                    tempSig.showTrough();
                elseif button=='s'
                    save([statsPath,'/stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                    disp('saving finished')
                elseif button=='r'
                    tempSig.showTrough();
                elseif button=='q'
                    return
                elseif button=='e'
                    %elseif button=='1'
                    tempSig.clearObj();
                    skipFlag=1;
                    break;
                elseif button=='c'
                    isClear=1;
                    save([statsPath,'/stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                    
                    disp('saving finished for previous droplets')
                end
                [x,y,button] = ginput(1);
            end
            if skipFlag
                skipFlag=0;
                ii=ii+1;
                disp(['skip_Trough: ', num2str(ii)]);
                %continue;
            end
        end
        
        %tempSig.showPeak();
        saveas(gcf,[imPath,'/raw/peak_',num2str(ii),'.png']);
        [dropTemp,osciTemp,rawTemp]=tempSig.output(ii);
        dropTemp.dropID(:)=tracks(ii).id;
        
        
        %{
        disp(length(tempSig.peakClass.time))
        if length(tempSig.peakClass.time)<peakThresh
            warning('number of peak too small, droplet discarded')
            ii=ii+1;
            tempSig.clearObj();
            continue;
        end
        %}
        rawTemp.time=frameStep*rawTemp.frame;

        eval(['rawTemp.int2=tracks(',num2str(ii),').Feat.',refNames{1},';']);
        
        osciTemp.int2=grpstats(rawTemp.int2(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
        %dropTemp.int2=mean(rawTemp.int2);
        
        %dropTemp.area=mean(areaSig);
        %dropTemp.meanMove=mean(sqrt(diff(xcoorSig).^2+diff(ycoorSig).^2));
        %dropTemp.frameStep=zeros(size(dropTemp.area))+frameStep;
        
        osciTemp.area=grpstats(rawTemp.areaSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
        osciTemp.xcoord=grpstats(rawTemp.xcoorSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
        osciTemp.ycoord=grpstats(rawTemp.ycoorSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
        osciTemp.frameStep=zeros(size(osciTemp.area))+frameStep;
        osciTemp.dropID(:)=tracks(ii).id;
        
        drop=[drop;dropTemp];
        osci=[osci;osciTemp];
        raw=[raw;rawTemp];
        tempSig.clearObj();
        % tempSig.lowPass();
        disp([num2str(ii),'/',num2str(size(tracks))])
        
        ii=ii+1;
    end
    drop.posIdx=zeros(size(drop.dropID))+posIdx;
    drop.position=zeros(size(drop.dropID))+position(posIdx);
    osci.posIdx=zeros(size(osci.dropID))+posIdx;
    osci.position=zeros(size(osci.dropID))+position(posIdx);
    %[osci.posIdx]=osciposIdx{:};
    %osci.posIdx=zeros(size(osci.area))+posIdx;
    raw.posIdx=zeros(size(raw.dropID))+posIdx;
    dropDat=[dropDat;drop];
    osciDat=[osciDat;osci];
    rawDat=[rawDat;raw];
    if isPlotStats
        plotStats([imPath,'/stats/'],osci,drop,frameStep);
    end
    
    posIdx=posIdx+1;
    save(fullfile(imPath,['stats_summary_',datestr(now,formatOut),'.mat']),'raw','osci','drop','-v7.3')
    drop=[];osci=[];raw=[];
    %posIdx=posIdx+1;
end
%% peak correction

dataPath='PATH_TO_RESULT_DIRECTORY\peakSelection';
savePath='PATH_TO_RESULT_DIRECTORY'; %save fig path
poi=[1];
correctOscillationData(dataPath, savePath, genPos, poi, sumTrackTable)