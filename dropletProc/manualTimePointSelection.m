%% load dataset
%load('id20211220_MOs_segment_20211224T102146.mat')
%load('id20211220_MOs_cleaning_20211224T112107.mat')

% use this script after data cleaning

%% params
% save dir of analysis results
dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20220422_Chk1';
%dataPath='/Volumes/qiongy-data/users/Gembu/results/20200923_Molpholino_Tune';
% name of sub-directory
plotType={'timeAlignment_start'};
% position for analysis
position=[16];   % all position to be processedposition=genPos;   % all position to be processed
posIdx=1;

% channel for peak detection
colNames={'CFPstdInt'};
sigScale=1;

mannualEdit=1;          % whether do mannual edition
oldF=0;                 % old naming format
loadSaveFile=0;         % load saved file?(same as dataPath)

% data save option
formatOut=30;
isClear=0;          % flag for clear and position

tempSig=manualSelectionClass();   
skipFlag=0;

for n=1:length(genPos)
    if position==genPos(n)
        savePath=fullfile(dataPath, plotType{1}, posList{1}{n});
    end    
end

%% load saved file
if loadSaveFile
    saveFileName=ls([savePath,'\selected_time_*.mat']);
    %saveFileName=ls([dataPath,'Stat_Pos',num2str(position(1)),'_*.mat']);
    if isempty (saveFileName)
        warning('no save file found')
        return;
    else
        saveFileName=saveFileName(end,:);
        disp(['loading ',saveFileName])
        load(fullfile(savePath,saveFileName),'osci','drop','raw')
        saveFile=load(fullfile(savePath,saveFileName),'tempSig','ii','posIdx','isClear');
        clearFlag=saveFile.isClear;
    end
end

%% main loop by the position
while posIdx <= length(position)
    selectedTP=[];
    selectedDROP=[];
    % overwrite position index if loadSaveFile was true
    if posIdx==1 && loadSaveFile==1
        posIdx=saveFile.posIdx;
    end
    
    
    % load data
    for n=1:size(posList{1},2)
        if strcmp(['Pos',num2str(position(posIdx))],posList{1}{n})
            relPos=n;
            savePath=fullfile(dataPath, plotType{1}, posList{1}{n});
            mkdir(savePath);
            disp(relPos);
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
            
            eval(['doi=tracks(',num2str(ii),').Feat.',colNames{1},';']); % data of interest
            
            areaSig=tracks(ii).Feat.area;
            xcoorSig=tracks(ii).Feat.xcoord;
            ycoorSig=tracks(ii).Feat.ycoord;
            tSig=tracks(ii).Feat.t;
            figure(1);
            plot(tSig,doi);
            hold on;
            tempSig.setOri(doi); % target signal. original and processed signal are preped. find peak trough with raw data.
%             tempSig.feature=table(tSig,madSig,areaSig,xcoorSig,ycoorSig);
%             tempSig.centerSig();
%             tempSig.fillOut(); % filloutliers with smooth thresh; doesn't have much effect
%             tempSig.genEnv(); 
%             tempSig.genPeak();
%             tempSig.attachPeak();
            %disp(length(tempSig.peakClass.time));
            
        end
        %% selection timepoint
        if mannualEdit==1
            [x,y,button] = ginput(1);
            tmpTP=round(x);
            tmpV=doi(tmpTP);
            scatter(tmpTP, tmpV,60,'bv');
            while button~='2'
                if button==1 %left click
                    %tempSig.userTPSelect(x);
                    selectedTP=[selectedTP, tmpTP];
                    selectedDROP=[selectedDROP, ii];
                elseif button==3 %right click
                    %tempSig.userDelPeak(x);
                    selectedTP=selectedTP(1:end-1);
                    selectedDROP=selectedDROP(1:end-1);
                elseif button=='s'
                    save([savePath,'/selected_time_',num2str(position(posIdx)),'.mat'],...
                        'selectedTP','selectedDROP','tempSig','ii','posIdx','isClear','-v7.3');
                    disp('saving finished');
                elseif button=='r'
                    % does this work?
                    tempSig.showPeak();
                elseif button=='q' %exit
                    return
                elseif button=='e' %skip
                    tempSig.clearObj();
                    skipFlag=1;
                    break;
                elseif button=='c'
                    % does this work?
                    isClear=1;
                    save([savePath,'/selected_time_',num2str(ii),'.mat'],...
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
        clf;
        tempSig.clearObj();
        % tempSig.lowPass();
        disp([num2str(ii),'/',num2str(size(tracks))])
        
        ii=ii+1;
    end
%     drop.posIdx=zeros(size(drop.dropID))+posIdx;
%     drop.position=zeros(size(drop.dropID))+position(posIdx);
%     osci.posIdx=zeros(size(osci.dropID))+posIdx;
%     osci.position=zeros(size(osci.dropID))+position(posIdx);
%     %[osci.posIdx]=osciposIdx{:};
%     %osci.posIdx=zeros(size(osci.area))+posIdx;
%     raw.posIdx=zeros(size(raw.dropID))+posIdx;
%     dropDat=[dropDat;drop];
%     osciDat=[osciDat;osci];
%     rawDat=[rawDat;raw];
%     if isPlotStats
%         plotStats([savePath,'/stats/'],osci,drop,frameStep);
%     end
    
    posIdx=posIdx+1;
    %save(fullfile(savePath,['stats_summary_',datestr(now,formatOut),'.mat']),'raw','osci','drop','-v7.3')
    drop=[];osci=[];raw=[];
    %posIdx=posIdx+1;
end