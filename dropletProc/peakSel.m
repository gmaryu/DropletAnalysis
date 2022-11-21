% for peak detection, trimming and mannual editing after preprocessing
% mouse left: add point
%       right: delpoint
%       middle: one step backward
%keyboard q: quit(without saving)
%         s: saving
%         c: saving all previous points
%         e: skip current point
%         r: refresh plot
dataPath='T:\Shiyuan\2019\11 05 2019 epi 4x\merged\';	%path of Pos<position>_bgCor_<time>.mat file from preprocess
position=[8];%all position to be processed
Channel2=1;% both YFP and mcherry present?
mannualEdit=1;%whether do mannual edition
oldF=0;% old naming format
loadSaveFile=1;%load saved file?(same as dataPath)
posIdx=1;

isSaveTrackFig=0;%save editing data to fig?
isSaveTrackPng=1;%save editing data to png?
isPlotStats=1;%plot statistics figure
frameStep=3;%min,frame rate

ampThresh=[-0.5,1];%binding range from smooth data to real data
adaptWindSize=2;%binding windoe size from smoothed data to real data

timeThresh=50;%Minimum time point for a track to be eligible for analysis

senThresh=0.05;%min prom of envolope difference
promSens=0.5;%added to envolope difference

smoothThresh=20;%window size for outlier detection

peakMinDist=30;%minimum peak distance
peakMaxDist=200; %/2 for generation of envolope
isClear=0;% flag for clear and position
formatOut=30;
filtThresh=[5,1,2];%max linear fitting(time vs amp) resid; max left amp/right amp; max prom/(max previous prom);
%robust linear fitting residual; left right amplitude ratio; max before/min
%amplitude
lbank=0.5/peakMinDist;
lpFilt = designfilt('lowpassiir','PassbandFrequency',lbank, ...
    'StopbandFrequency',lbank+0.05,'PassbandRipple',0.5, ...
    'StopbandAttenuation',65);%low pass filter reserved
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
if loadSaveFile
    saveFileName=ls([dataPath,'Stat_Pos',num2str(position(1)),'_*.mat']);
    if isempty (saveFileName)
        warning('no save file found')
        return;
    else
        saveFileName=saveFileName(end,:);
        disp(['loading ',saveFileName])
        load([dataPath,saveFileName],'osci','drop','raw')
        saveFile=load([dataPath,saveFileName],'tempSig','ii','posIdx','isClear');
        clearFlag=saveFile.isClear;
    end
end
while posIdx<=length(position)
    
    if posIdx==1&&loadSaveFile==1
        posIdx=saveFile.posIdx;
        
    end
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
    imPath=[dataPath,'Pos',num2str(position(posIdx)),'_figure'];
    mkdir(imPath)
    mkdir([imPath,'\raw'])
    mkdir([imPath,'\stats'])
    if size(fileName,1)~=1
        error('Check the number of file specified')
    end
    load([dataPath,fileName],'resultTable')
    %cond='wee1=1';
    %resultTable(resultTable.t>500,:)=[];

    resultTable=resultTable(resultTable.track>=1,:);
    idlist=unique(resultTable.track);
    ii=1;
    
    while ii<=length(unique(resultTable.track))
        if ii==1&&loadSaveFile==1&&~clearFlag
            ii=saveFile.ii;
            
            tempSig=saveFile.tempSig;
            
        else
            if ii==1&&loadSaveFile==1
                ii=saveFile.ii;
            end
            dellist=[];
            if Channel2 &&~oldF
                madSig=resultTable.CFPstdInt(resultTable.track==idlist(ii));
                aveSig=resultTable.CFPdivInt(resultTable.track==idlist(ii));               
                
            elseif Channel2
                madSig=resultTable.C1madInt(resultTable.track==idlist(ii));
                aveSig=resultTable.C1divInt(resultTable.track==idlist(ii));               
               
            else
                madSig=resultTable.CFPstdInt(resultTable.track==idlist(ii));
                aveSig=resultTable.divInt(resultTable.track==idlist(ii));
            end
            areaSig=resultTable.area(resultTable.track==idlist(ii));
            xcoorSig=resultTable.xcoord(resultTable.track==idlist(ii));
            ycoorSig=resultTable.ycoord(resultTable.track==idlist(ii));
            tSig=resultTable.t(resultTable.track==idlist(ii));
            
            if length(tSig)<timeThresh
                ii=ii+1;
                tempSig.clearObj();
                continue;
            end
%             if mean(resultTable.mCherry(resultTable.track==idlist(ii)))<1500
%                 
%                 ii=ii+1;
%                 continue;
%             end
            if Channel2&&~oldF

            tempSig.setRef(resultTable.mCherry(resultTable.track==idlist(ii)),... 
                resultTable.YFP(resultTable.track==idlist(ii)));
            elseif Channel2
             divInt.setRef(resultTable.C1divInt(resultTable.track==idlist(ii)),...
                resultTable.C2divInt(resultTable.track==idlist(ii)));
               
            end
               % resultTable.DAPI(resultTable.track==idlist(ii)));
            tempSig.setOri(aveSig);
            tempSig.feature=table(tSig,madSig,areaSig,xcoorSig,ycoorSig);
            tempSig.centerSig();
            tempSig.fillOut();
            tempSig.genEnv();
            tempSig.genPeak();
            tempSig.attachPeak();
            if length(tempSig.peakClass.time)<peakThresh
                ii=ii+1;
                tempSig.clearObj();
                continue;
            end
            tempSig.genTrough()
            tempSig.checkPeak();
        end
        tempSig.showPeak();
        
        [x,y,button] = ginput(1);
        
        while button~='2'
            if button==1
                tempSig.userAddPeak(x);
            elseif button==3
                tempSig.userDelPeak(x);
            elseif button==2
                tempSig.peakClass=tempSig.backClass;
                tempSig.showPeak();
            elseif button=='s'
                save([dataPath,'Stat_',fileName(1:end-4),'.mat'],...
                    'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                disp('saving finished')
                
            elseif button=='r'
                tempSig.showPeak();
            elseif button=='q'
                return
            elseif button=='1'
                tempSig.clearObj();
                skipFlag=1;
                break;
            elseif button=='c'
                isClear=1;
                save([dataPath,'Stat_',fileName(1:end-4),'.mat'],...
                    'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                disp('saving finished for previous droplets')
                
            end
            [x,y,button] = ginput(1);
        end
        if skipFlag
            skipFlag=0;
            ii=ii+1;
            continue;
        end
        tempSig.genTrough();
        tempSig.showTrough();
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
                save([dataPath,'Stat_',fileName(1:end-4),'.mat'],...
                    'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                disp('saving finished')
            elseif button=='r'
                tempSig.showTrough();
            elseif button=='q'
                return
            elseif button=='1'
                tempSig.clearObj();
                skipFlag=1;
                break;
            elseif button=='c'
                isClear=1;
                save([dataPath,'Stat_',fileName(1:end-4),'.mat'],...
                    'osci','drop','raw','tempSig','ii','posIdx','isClear','-v7.3');
                disp('saving finished for previous droplets')
            end
            [x,y,button] = ginput(1);
        end
        if skipFlag
            skipFlag=0;
            ii=ii+1;
            continue;
        end
        
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
            saveas(gcf,[imPath,'\raw\peak_',num2str(idlist(ii)),'.fig'])
        end
        if isSaveTrackPng
            saveas(gcf,[imPath,'\raw\peak_',num2str(idlist(ii)),'.png'])
        end
        [dropTemp,osciTemp,rawTemp]=tempSig.output(ii);
        
        rawTemp.time=frameStep*rawTemp.frame;
        
        if Channel2&&~oldF
            rawTemp.int2=resultTable.mCherry(resultTable.track==idlist(ii));
            rawTemp.int3=resultTable.YFP(resultTable.track==idlist(ii));
            osciTemp.int2=grpstats(rawTemp.int2(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            osciTemp.int3=grpstats(rawTemp.int3(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            dropTemp.int3=mean(rawTemp.int3);            
            dropTemp.int2=mean(rawTemp.int2);
        elseif Channel2
            rawTemp.int2=resultTable.C2divInt(resultTable.track==idlist(ii));
            osciTemp.int2=grpstats(rawTemp.int2(rawTemp.peakCount>0),rawTemp.peakCount(rawTemp.peakCount>0));
            dropTemp.int2=mean(rawTemp.int2);
            
   
        end
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
        %        tempSig.lowPass();
        ii=ii+1;
        disp([num2str(ii),'/',num2str(length(unique(resultTable.track)))])
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
    save([dataPath,fileName(1:end-4),'_Stat.mat'],'raw','osci','drop','-v7.3')
    
    drop=[];osci=[];raw=[];
end

save([dataPath,'Stats_Pos',num2str(position),'_',datestr(now,formatOut),'.mat'])
