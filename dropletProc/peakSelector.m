% for peak detection, trimming and mannual editing after preprocessing
% mouse left: add point
%       right: delpoint
%       middle: one step backward
%keyboard q: quit(without saving)
%         s: saving
%         c: saving all previous points
%         e: skip current point
%         r: refresh plot

% use this script after data cleaning

%% params

% save dir of analysis results
dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20200224_mRNA_inhibitors';

% name of sub-directory
plotType={'peakSelection'};

% position for analysis
position=[19];%all position to be processed
%position=[2,3,5,7];%all position to be processed

% for i=1:length(genPos)
%     posList{1}{i}=['Pos',num2str(genPos(i))];
% end

% channel for peak detection
colNames={'FRET','FRET'};% both YFP and mcherry present?
refNames={'mCherry','CFPstdInt'};
%colName=1;% both YFP and mcherry present?
posIdx=1;

mannualEdit=1;%whether do mannual edition
oldF=0;% old naming format
loadSaveFile=1;%load saved file?(same as dataPath)

% data save option
isSaveTrackFig=0;%save editing data to fig?
isSaveTrackPng=1;%save editing data to png?
isPlotStats=0;%plot statistics figure
formatOut=30;

% params for peak detection
frameStep=3;%min,frame rate
ampThresh=[-0.5,1];%binding range from smooth data to real data
adaptWindSize=2;%binding windoe size from smoothed data to real data
%timeThresh=50;%Minimum time point for a track to be eligible for analysis
senThresh=0.05;%min prom of envolope difference
promSens=0.5;%added to envolope difference
smoothThresh=20;%window size for outlier detection
peakMinDist=30;%minimum peak distance
peakMaxDist=200; %/2 for generation of envolope

isClear=0;% flag for clear and position
filtThresh=[5,1,2];%max linear fitting(time vs amp) resid; max left amp/right amp; max prom/(max previous prom);

%robust linear fitting residual; left right amplitude ratio; max before/min
%amplitude
lbank=0.5/peakMinDist;
lpFilt = designfilt('lowpassiir','PassbandFrequency',lbank, ...
    'StopbandFrequency',lbank+0.05,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65);% low pass filter reserved
% baseFilt = designfilt('lowpassiir','PassbandFrequency',0.005, ...
%     'StopbandFrequency',0.025,'PassbandRipple',1, ...
%     'StopbandAttenuation',10);
tempSig=signalClass(senThresh,promSens,smoothThresh,peakMinDist,filtThresh,...
    lpFilt,peakMaxDist,adaptWindSize,ampThresh);%temp droplet class
trackCell=[];
drop=[];
osci=[];
raw=[];
peakThresh=3;
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
try
    while posIdx <= length(position)
        % overwrite position index if loadSaveFile was true
        if posIdx==1 && loadSaveFile==1
            posIdx=saveFile.posIdx;
        end
        
        % load data
        for n=1:size(posList{1},2)
            if strcmp(['Pos',num2str(position(posIdx))],posList{1}{n})
                relPos=n;
                %disp(relPos);
                % for save results
                imPath=fullfile(dataPath, plotType{1}, posList{1}{n});
                mkdir(imPath)
                rawPath=fullfile(imPath,'\raw');
                mkdir(rawPath);
                statsPath=fullfile(imPath,'\stats');
                mkdir(statsPath);
                break;
            end
        end
        
        %{
    fileName=ls([dataPath,'Pos',num2str(position(posIdx)),'_bgCor_*']);
    if size(fileName,1)==2
        fileList=contains(string(fileName),'_Stat');
        fileName=fileName(~fileList,:);
        while ~isempty(fileName)
            if strcmp(fileName(end),' ')
                fileName(end)=[];
            else
                break;
            end
        end
    end
    
    if size(fileName,1)~=1
        error('Check the number of file specified')
    end
        %}
        
        tracks=sumTrackTable{n};
        
        %load([dataPath,fileName],'resultTable')
        %cond='wee1=1';
        %resultTable(resultTable.t>500,:)=[];
        
        %resultTable=resultTable(resultTable.track>=1,:);
        %idlist=unique(resultTable.track);
        ii=1;
        
        while ii <= length(tracks)
            disp(['ii: ',num2str(ii)]);
            if ii==1&&loadSaveFile==1&&~clearFlag
                disp(saveFile.ii);
                ii=saveFile.ii;
                tempSig=saveFile.tempSig;
                %disp(['1,load,clear ',num2str(ii)]);
            else
                if ii==1 && loadSaveFile==1
                    ii=saveFile.ii;
                    %disp(ii);
                    %disp(['1,load  ',num2str(ii)]);
                end
                eval(['aveSig=tracks(',num2str(ii),').Feat.',colNames{1},';']);
                eval(['madSig=tracks(',num2str(ii),').Feat.',colNames{2},';']);
                
                %{
%             dellist=[];
%             if colName &&~oldF
%                 madSig=resultTable.CFPstdInt(resultTable.track==idlist(ii));
%                 aveSig=resultTable.CFPdivInt(resultTable.track==idlist(ii));
%
%             elseif colName
%                 madSig=resultTable.C1madInt(resultTable.track==idlist(ii));
%                 aveSig=resultTable.C1divInt(resultTable.track==idlist(ii));
%
%             else
%                 madSig=resultTable.CFPstdInt(resultTable.track==idlist(ii));
%                 aveSig=resultTable.divInt(resultTable.track==idlist(ii));
%             end
                %}
                areaSig=tracks(ii).Feat.area;
                xcoorSig=tracks(ii).Feat.xcoord;
                ycoorSig=tracks(ii).Feat.ycoord;
                tSig=tracks(ii).Feat.t;
                %{
            % thresholding by timeThresh has done
%             if length(tSig)<timeThresh
%                 ii=ii+1;
%                 tempSig.clearObj();
%                 continue;
%             end
            
%             if mean(resultTable.mCherry(resultTable.track==idlist(ii)))<1500
%
%                 ii=ii+1;
%                 continue;
%             end
                %}
                tempSig.setRef(eval(['tracks(',num2str(ii),').Feat.',refNames{1},';']),...
                    eval(['tracks(',num2str(ii),').Feat.',refNames{2},';']));
                tempSig.setOri(aveSig);
                tempSig.feature=table(tSig,madSig,areaSig,xcoorSig,ycoorSig);
                tempSig.centerSig();
                tempSig.fillOut();
                tempSig.genEnv();
                tempSig.genPeak();
                tempSig.attachPeak();
                disp(length(tempSig.peakClass.time));
                if length(tempSig.peakClass.time)<peakThresh
                    ii=ii+1;
                    tempSig.clearObj();
                    disp('below peakThresh');
                    continue;
                end
                
                tempSig.genTrough();
                tempSig.checkPeak();
            end
            tempSig.showPeak();
            
            %% selection peaks
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
                    save([statsPath,'\stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
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
                    save([statsPath,'\stats_track_',num2str(ii),'.mat'],...
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
            
            tempSig.genTrough();
            tempSig.showTrough();
            
            %% trough selection
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
                    save([statsPath,'\stats_track_',num2str(ii),'.mat'],...
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
                    save([statsPath,'\stats_track_',num2str(ii),'.mat'],...
                        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                    %}
                    disp('saving finished for previous droplets')
                end
                [x,y,button] = ginput(1);
            end
            if skipFlag
                skipFlag=0;
                ii=ii+1;
                disp(['skip_Trough: ', num2str(ii)]);
                continue;
            end
            disp('move to following process');
            %%
            if length(tempSig.troughClass.time)~=length(tempSig.peakClass.time)-1
                tempSig.genTrough();
                warning('trough detection number do not match, manual data discarded!')
            end
            if length(tempSig.peakClass.time)<peakThresh
                warning('number of peak too small, droplet discarded')
                ii=ii+1;
                tempSig.clearObj();
                continue;
            end
            if isSaveTrackFig
                saveas(gcf,[imPath,'\raw\peak_',num2str(ii),'.fig'])
            end
            if isSaveTrackPng
                saveas(gcf,[imPath,'\raw\peak_',num2str(ii),'.png'])
            end
            [dropTemp,osciTemp,rawTemp]=tempSig.output(ii);
            
            rawTemp.time=frameStep*rawTemp.frame;
            
            %{
        if colNames&&~oldF
            %rawTemp.int2=resultTable.mCherry(resultTable.track==idlist(ii));
            %rawTemp.int3=resultTable.YFP(resultTable.track==idlist(ii));
            eval(['rawTemp.int2=tracks(',num2str(n),').Feat.',refNames{1},';']);
            eval(['rawTemp.int2=tracks(',num2str(n),').Feat.',refNames{2},';']);
            osciTemp.int2=grpstats(rawTemp.int2(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            osciTemp.int3=grpstats(rawTemp.int3(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            dropTemp.int3=mean(rawTemp.int3);
            dropTemp.int2=mean(rawTemp.int2);
        elseif colNames
            %}
            eval(['rawTemp.int2=tracks(',num2str(ii),').Feat.',refNames{1},';']);
            %rawTemp.int2=resultTable.C2divInt(resultTable.track==idlist(ii));
            osciTemp.int2=grpstats(rawTemp.int2(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            dropTemp.int2=mean(rawTemp.int2);
            
            dropTemp.area=mean(areaSig);
            dropTemp.meanMove=mean(sqrt(diff(xcoorSig).^2+diff(ycoorSig).^2));
            dropTemp.frameStep=zeros(size(dropTemp.area))+frameStep;
            
            osciTemp.area=grpstats(rawTemp.areaSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            osciTemp.xcoord=grpstats(rawTemp.xcoorSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            osciTemp.ycoord=grpstats(rawTemp.ycoorSig(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            osciTemp.frameStep=zeros(size(osciTemp.area))+frameStep;
            
            drop=[drop;dropTemp];
            osci=[osci;osciTemp];
            raw=[raw;rawTemp];
            tempSig.clearObj();
            % tempSig.lowPass();
            disp([num2str(ii),'/',num2str(size(tracks))])
            ii=ii+1;
        end
        
        drop.posIdx=zeros(size(drop.area))+posIdx;
        osci.posIdx=zeros(size(osci.area))+posIdx;
        raw.posIdx=zeros(size(raw.time))+posIdx;
        dropDat=[dropDat;drop];
        osciDat=[osciDat;osci];
        rawDat=[rawDat;raw];
        if isPlotStats
            plotStats([imPath,'\stats\'],osci,drop,frameStep);
        end
        
        posIdx=posIdx+1;
        save(fullfile(imPath,['stats_summary_',datestr(now,formatOut),'.mat']),'raw','osci','drop','-v7.3')
        
        drop=[];osci=[];raw=[];
    end
catch
    isClear=1;
    save([statsPath,'\stats_track_',num2str(ii),'.mat'],...
        'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
    disp('bad mouse control')
end