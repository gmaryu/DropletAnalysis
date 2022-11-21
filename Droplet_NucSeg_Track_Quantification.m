clear all

%% global params
addpath(pwd)
dateList={'DIRECTORY_NAME'};
ds = split(dateList, '_');
datadate=['id',ds{1}];
genPos=[0,1,2,3,4,5,6,7];
for i=1:length(genPos)
    posList{1}{i}=['Pos',num2str(genPos(i))];
end

frameNum=300; % number of frames for analysis (this params also works as list with posList)
extraFrameNum=0; % additional analysis as a snapshot?

resAmp=1.0; % resolution
minDopDia=30*resAmp; % minimum drop diameter?
% backThresh=500; % background threshold?
EccentricityT=0.95;
prthresh=0.35; % perimeter thresh?
lenThresh=350; % minimum tracking length
areaT=[1e3*resAmp^2, 1e4*resAmp^2]; % min and max threshold of area size
areaChangePercent=0.5; % how much area change is allowed
formatOut=30; % date format for 'yyyymmddTHHMMSS'
%   nucThresh=1;

nameFormat='img_000000%03d_4-BF_000.tif';
cFormat{1}='img_000000%03d_5-CFP_000.tif';
cFormat{2}='img_000000%03d_8-Custom_000.tif';
cFormat{3}='img_000000%03d_1-DAPI_000.tif';
cName={'CFP','FRET','BFP'};

%% segmentation and intensity quantification
for pp=1:length(dateList)
    nf=frameNum(pp);
    sumSegTable = cell(1,length(genPos));
    sumIntTable_n = cell(1,length(genPos));
    sumIntTable_c = cell(1,length(genPos));
    sunIntTable_w = cell(1,length(genPos));
    sumExtTable = cell(1,length(genPos));
    sumResultTable = cell(1,length(genPos));
    dataPathTable = cell(1,length(genPos));
    %for pj=1:length(posList{pp})
    parfor pj=1:length(posList{pp})
        %disp(genPos(pj));
        %%
        foldername=['DATA_PARENT_PATH',dateList{pp},'/',posList{pp}{pj}];
        
        addpath(genpath(pwd))
        cd(foldername)
        
        resultTable=[];
        exTable=[];
        segChTable = cell(1, nf);
        intChTable_n = cell(1, nf);
        intChTable_c = cell(1, nf);
        intChTable_w = cell(1, nf);
        areaTable = cell(1, nf);
        extChTable = cell(1, nf);
        %for i=1:25:nf
        for i=1:nf
            %i=125;
            X = ['Pos',num2str(genPos(pj)), ': ', num2str(i),' Frame'];
            %disp(X);
            
            % load bright field image
            fname= sprintf(nameFormat,i-1);
            X = [X, ' / ', fname];
            disp(X);
            img= imread(fname);
            %% segmentation with BF image
            [sg, S1]=BFsegmentationDroplets(img,i,resAmp,areaT,EccentricityT,prthresh,cFormat);
            
            % add skewness and stdInt to Structure
            S1=skwDroplet(S1);
            
            %%
            %{
            % test code for label and regionprops is coupled.
            smplist=randi([1,size(S1,1)],1,10);
            disp(smplist);
            for n=1:length(smplist)
                tar=S1(smplist(n)).Label;
                cnt=S1(smplist(n)).Centroid;
                disp(cnt);
                figure; imagesc(sg==tar);
                hold on;
                scatter(cnt(1), cnt(2), 150,'d');
            end
            %}
            
            %% detection of nuclear ---------------
            fname_mk= sprintf(cFormat{1},i);
            img2=double(imread(fname_mk)); % nuc marker image
            
            %% detection of nuclear ---------------
            fname_mk= sprintf(cFormat{1},i-1);
            img2=double(imread(fname_mk)); % nuc marker image
            
            %-- params for nuc segmentation --
            cvth1=0.13;
            cvth2=[1.2, 1.0, 0.80, 0.4, 0.25, 0.14];
            stdth=100;
            skwth=0;
            % --
            
            [imnc,imcy,imnc2,imcy2]=nucDetector(sg,img2,S1,cvth1,cvth2,stdth,skwth);
            
            % save nuc-image
            savedir=['PATH_TO_RESULT_DIRECTORY',dateList{pp},'/nuc_cyt_seg/',posList{pp}{pj},'/'];
            nameFormat1='img_000000%03d_nuc_000.jpg';
            fname1= sprintf(nameFormat1,i-1);
            nameFormat1='img_000000%03d_cyto_000.jpg';
            fname2= sprintf(nameFormat1,i-1);
            if ~exist(savedir, 'dir')
                mkdir(savedir)
            end
            imwrite(imnc2, [savedir,fname1]);
            imwrite(imcy2, [savedir,fname2]);            
            
            %% quantify for resultmat
            S2=regionprops(sg,img2,'MeanIntensity','Centroid','Area',...
                'Perimeter','BoundingBox','PixelIdxList','PixelValues','Eccentricity','EulerNumber');
            delid2=[S2.Area]==0;
            S2(delid2==1)=[];
            
            cent=reshape([S2.Centroid],2,length([S2.Area]));
            T=table(zeros(size([S2.MeanIntensity]'))+i,cent(1,:)',cent(2,:)',...
                [S2.Area]',[S2.Perimeter]',[S1.stdInt]',[S2.MeanIntensity]', [S1.skewness]');
            T.Properties.VariableNames = {'t','xcoord', 'ycoord', 'area', 'Perimeter',...
                [cName{1},'stdInt'],[cName{1},'meanInt'], 'skewness'};
            segChTable{i}=T;
            
           
           % quantify nuclear area size
            S3=regionprops(imnc2,img2,'Area');
            S3(delid2==1)=[];
            S4=regionprops(imcy2,img2,'Area');
            S4(delid2==1)=[];
            T2=table([S3.Area]',[S4.Area]');
            T2.Properties.VariableNames = {'Nuc_area', 'Cyto_area'};
            areaTable{i}=T2;
            
            
            %{
            if i<extraFrameNum
                %exT=[T,ttemp];
                exT=[T,exTab];
                exTable=[exTable;exT];
            end
            %resultTable=[resultTable;T(2:end,:)];
            %}
            %% record other channel's intensity
            intTab_n= cell(1,length(cFormat));
            intTab_c= cell(1,length(cFormat));
            intTab_w= cell(1,length(cFormat));

            for ijk=1:length(cFormat)
                fname= sprintf(cFormat{ijk},i-1);
                C=imread(fname);
                C=uint16(C);
                C=double(C);
                S2n=regionprops(imnc2,C,'MeanIntensity','PixelIdxList','PixelValues','Area');
                S2n(delid2==1)=[];
                S2c=regionprops(imcy2,C,'MeanIntensity','PixelIdxList','PixelValues','Area');
                S2c(delid2==1)=[];
                S2w=regionprops(sg,C,'MeanIntensity','PixelIdxList','PixelValues','Area');
                S2w(delid2==1)=[];
                ttemp_n=table([S2n.MeanIntensity]'); % intensity of other channel
                ttemp_c=table([S2c.MeanIntensity]'); % intensity of other channel
                ttemp_w=table([S2w.MeanIntensity]'); % intensity of other channel
                %disp(size(ttemp));
                ttemp_n.Properties.VariableNames={[char(cName(ijk)),'_nuc']};
                ttemp_c.Properties.VariableNames={[char(cName(ijk)),'_cyto']};
                ttemp_w.Properties.VariableNames={[char(cName(ijk)),'_whole']};
                %sumResultTable{pj,2,nf}=ttemp;
                intTab_n{ijk}=ttemp_n;
                intTab_c{ijk}=ttemp_c;
                intTab_w{ijk}=ttemp_w;
            end
            intChTable_n{i}=intTab_n;
            intChTable_c{i}=intTab_c;
            intChTable_w{i}=intTab_w;
            %disp(size(ttemp));
            
            if i<extraFrameNum
                fname= sprintf(cFormat{end},i-1);
                C=imread(fname);
                C=uint16(C);
                C=double(C);
                S11=regionprops(sg,C,'MeanIntensity','PixelIdxList','PixelValues');
                S11(delid==1)=[];
                ttemp=table([S11.MeanIntensity]');
                ttemp.Properties.VariableNames=cName(end);
                extChTable{i}=ttemp;
            
            end
            %}
        end
        sumIntTable_n{pj}=intChTable_n;
        sumIntTable_c{pj}=intChTable_c;
        sumIntTable_w{pj}=intChTable_w;
        sumAreaTable{pj}=areaTable;
        sumSegTable{pj}=segChTable;
        sumExtTable{pj}=extChTable;
        dataPathTable{pj}=foldername;
    end
    
    %% gather quantified data for multiple position into single varialbe
    for pj =1:length(posList{pp}) %position
        
        tmpSegTable=sumSegTable{pj}; % 1st position
        tmpIntTable_n=sumIntTable_n{pj};
        tmpIntTable_c=sumIntTable_c{pj};
        tmpIntTable_w=sumIntTable_w{pj};
        tmpAreaTable=sumAreaTable{pj};
        
        for i=1:length(tmpSegTable) % frame number
            % summlize other channels data to one table
            for k=1:length(tmpIntTable_n{i}) % channel
                if k==1
                    allIntTable= tmpIntTable_n{i}{1};
                    allIntTable= horzcat(allIntTable, tmpIntTable_c{i}{1});
                    allIntTable= horzcat(allIntTable, tmpIntTable_w{i}{1});
                else
                    allIntTable= horzcat(allIntTable, tmpIntTable_n{i}{k});
                    allIntTable= horzcat(allIntTable, tmpIntTable_c{i}{k});
                    allIntTable= horzcat(allIntTable, tmpIntTable_w{i}{k});
                end
            end
            allIntTable= horzcat(allIntTable, tmpAreaTable{i});
            tmpIntTable{i}= allIntTable;
            if i == 1
                posResultTable=tmpSegTable{i};
                posIntTable2=tmpIntTable{i};
            else
                posResultTable=vertcat(posResultTable, tmpSegTable{i});
                posIntTable2=vertcat(posIntTable2, tmpIntTable{i});
            end
        end
        sumResultTable{pj}=horzcat(posResultTable, posIntTable2);
    end
    
end
save([datadate,'_NucSeg_segment_',datestr(now,formatOut)],'-v7.3')

%% tracking
sumTrackTable=cell(1,length(genPos));
for pos=1:length(sumResultTable)
    %parfor pos=1:length(sumResultTable)
    resultTable=sumResultTable{pos};
    % refernce channel
    refCh=cName{1};
    refInt=eval(['resultTable.',refCh,'meanInt']);
    % tracking params
    par=zeros(14,1);
    par(1)=180*resAmp*0.5; % max dist 60
    par(2)=60*resAmp; % std dist 15
    par(3)=std(log10(resultTable.area)); % log scale area
    par(4)=min(mad(refInt)/3, 10000);
    %par(4)=min(mad(resultTable.mCherrymeanInt)/3, 10000);% mad-> mean or median OR constant value
    %par(5)=min(mad(resultTable.YFP)/1.5,5);
    par(5)=999999;
    par(6:9)=1;
    par(10:13)=0;
    
    % tracking
    %[tracks,trackArray]=CellTrackerNostat(resultTable,par);
    [tracks,trackArray]=CellTrackerNostat(resultTable,refCh,par);
    delArray=trackArray;
    filtlab1=[tracks.Len]<lenThresh; % automatic filtering by label
    ar=zeros(size(filtlab1));
    for i=1:length(tracks)
        ar(i)=range([tracks(i).Feat.area])/mean([tracks(i).Feat.area]);
        if ar(i)>areaChangePercent
            filtlab1(i)=1; % add to delete label vector
        end
    end
    
    for i=1:max(trackArray)
        if filtlab1(i)
            delArray(trackArray==i)=0;
        end
    end
    resultTable.track=trackArray;
    resultTable.delArray=delArray;
    sumResultTable{pos}=resultTable;
    sumTrackTable{pos}=tracks;
    %sumTrackTable{pos}=resultTable;
end
save([datadate,'_NucSeg_trackinig_',datestr(now,formatOut)],'sumResultTable','sumTrackTable')

%% data cleaning
for i=1:size(sumTrackTable,2)
    tracks=sumTrackTable{i};
    %% cut short trace
    tracks=cutShortTrace(tracks,nf,0.9);

    %% Calc ratiometric value
    tracks=calcRatio(tracks,'FRET_whole','CFP_whole','FC_ratio_whole');

    %% with or without Nuc classification
    tracks=classifyNucCyt(tracks);
    sumTrackTable{i}=tracks;
end
save([datadate,'_NucSeg_cleaning_',datestr(now,formatOut)],'sumResultTable','sumTrackTable','dateList','posList','genPos')
%mvSegmentationResult(dateList, posList)

%% plot result 
plotParam{1}=sumTrackTable;
plotParam{2}=posList;
plotParam{3}=dateList;
plotParam{4}=genPos;
plotParam{5}='PATH_TO_RESULT_DIRECTORY'; % store directory
plotParam{6}={'FC_ratio_std'}; %plotType
% left Y-axis
plotParam{7}={'FC_ratio_whole'}; % colum name for left y-axis
plotParam{8}='FRET/CFP Ratio'; % label for left y-axis
% right Y-axis 
plotParam{9}={'CFPstdInt'}; % colum name for right y-axis
plotParam{10}='CFP intensity Std';% label for right y-axis
% time
plotParam{11}=4; % interval
plotParam{12}='Time (min)'; % label for x-axis
plotParam{13}='all';
plotParam{14}='png';
plotColData(plotParam);
