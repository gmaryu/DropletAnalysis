clear all
%{

%}
%% global params
addpath(pwd)
dateList={'DIRECTORY_NAME'}; % YOUR DATA PATH
ds = split(dateList, '_');
datadate=['id',ds{1}];
genPos=[0,1,2,3,4,5,6,7]; % YOUR DATA POSITION
for i=1:length(genPos)
    posList{1}{i}=['Pos',num2str(genPos(i))];
end

frameNum=200; % NUMBER OF FRAMES FOR ANALYSIS (this params also works as list with posList)
extraFrameNum=0; % additional analysis as a snapshot?

resAmp=1.0; % resolution 1x1-> 1.0; 2x2->0.5
minDopDia=30*resAmp; % minimum drop diameter?
% backThresh=500; % background threshold?
EccentricityT=0.95;
prthresh=0.35; % perimeter thresh?
lenThresh=350; % minimum tracking length
areaT=[1e3*resAmp^2, 1e4*resAmp^2]; % min and max threshold of area size
areaChangePercent=0.5; % how much area change is allowed
formatOut=30; % date format for 'yyyymmddTHHMMSS'
%   nucThresh=1;

% IMAGE FILE NAME FOR QUANTIFICATION
nameFormat='img_000000%03d_4-BF_000.tif';
cFormat{1}='img_000000%03d_5-CFP_000.tif';
cFormat{2}='img_000000%03d_8-Custom_000.tif';
cFormat{3}='img_000000%03d_6-YFP_000.tif';
cName={'CFP','FRET','YFP'};

%% segmentation and intensity quantification
for pp=1:length(dateList)
    nf=frameNum(pp);
    sumSegTable = cell(1,length(genPos));
    sunIntTable_w = cell(1,length(genPos));
    sumExtTable = cell(1,length(genPos));
    sumResultTable = cell(1,length(genPos));
    dataPathTable = cell(1,length(genPos));
    %for pj=1:length(posList{pp})
    parfor pj=1:length(posList{pp})
        %disp(genPos(pj));
        foldername=['DATA_PARENT_PATH',dateList{pp},'/',posList{pp}{pj}];
        addpath(genpath(pwd))
        cd(foldername)
        
        resultTable=[];
        exTable=[];
        segChTable = cell(1, nf);
        intChTable_w = cell(1, nf);
        areaTable = cell(1, nf);
        extChTable = cell(1, nf);
        %for i=1:25:nf
        for i=1:nf
            X = ['Pos',num2str(genPos(pj)), ': ', num2str(i),' Frame'];
            %disp(X);
            %% segmentation with BF image
            % load bright field image
            fname= sprintf(nameFormat,i-1);
            X = [X, ' / ', fname];
            disp(X);
            img= imread(fname);
            gimg= double(img);
            
            % normalize
            ngimg=(gimg-min(gimg(:)))/(max(gimg(:))-min(gimg(:)));
            
            % edge detection for segmentation
            c=0.0005;
            gimg2=imgaussfilt(ngimg,6)-imgaussfilt(ngimg,2);
            gimg3=gimg2;
            gimg3(ngimg<0)=0;
            
            sigma=3;
            [gDxx,gDxy,gDyy] = Hessian2D(gimg2,sigma);
            % Correct for scale
            gDxx = (sigma^2)*gDxx;
            gDxy = (sigma^2)*gDxy;
            gDyy = (sigma^2)*gDyy;
            
            % Calculate (abs sorted) eigenvalues and vectors
            [gLambda2,gLambda1,~,~]=eig2image(gDxx,gDxy,gDyy);
            gS2 = gLambda1.^2;% + Lambda2.^2;
            gI2=(ones(size(gimg2))-exp(-gS2/c)).*(gLambda1<0).*gimg3;
            tmpthresh=graythresh(gI2)/4;
            ngimg=gI2>tmpthresh;
            ngimgl=bwlabel(ngimg); %binalize and labling
            %figure; imshow(ngimgl);
            %saveas(gcf,'ngimgl.tif');
            
            gS=regionprops(ngimgl,'Area','Eccentricity');
            delidx=[gS.Area]<(50*resAmp^2)|(([gS.Area]<1000*resAmp^2)&[gS.Eccentricity]<0.95); % pick-up of junk labels
            
            % replace labels to zero
            for iij=1:length(delidx)
                if delidx(iij)
                    ngimgl(ngimgl==iij)=0;
                end
            end
            %figure; imshow(ngimgl);
            ngimg=logical(ngimgl); % binalize of label image
            gimg2=imclose(ngimg,strel('disk',round(5*resAmp))); % morphology filter for binalized image
            gimgo=bwareaopen(gimg2,1500); % delete small objects
            %{bainari dakede segmentation dekirunodeha?}
            ngimg=gimg3&ngimg;
            seeds=imextendedmax(bwdist(ngimg),2);
            gimg4=max(gimg(:))- gimg;
            gss1=watershed(imimposemin(gimg4,seeds));
            gss2=watershed(imimposemin(gimg4,gss1==0|seeds));
            %figure; imagesc(gss2);
            gss2(gss2==1)=0;
            gss3=logical(gss2);
            gss3=imopen(gss3,strel('disk',round(2*resAmp)));
            gss2=double(gss2).*double(gss3);
            
            fname= sprintf(cFormat{1},i-1);
            B=imread(fname);
            gB=double(B);
            S1=regionprops(gss2,B,'Centroid','Area','Perimeter',...
                'PixelIdxList','PixelValues','Eccentricity','EulerNumber');
            
            % collection and overwrite of junk label
            %disp(size(S1));
            delid=[S1.Area]>areaT(2)|[S1.Area]<areaT(1); % too big or too small
            delid=delid|[S1.Eccentricity]>EccentricityT; % roundness
            delid=delid|[S1.Perimeter].^2./[S1.Area]>(1+prthresh)*4*pi; % mismutch of perimeter and area ?
            delid=delid|[S1.EulerNumber]>1|[S1.EulerNumber]<1;
            %disp(sum(delid));
            S1(delid==1)=[];
            
            % apply delid to watershed result (cleaned segmentation result)
            %figure; imagesc(gss2);
            ss5=gss2;
            for ij=1:max(ss5(:))
                %for ij=2:max(gss2(:))
                if delid(ij)
                    ss5(ss5==ij)=0;
                end
            end
            
            ss5l=unique(ss5);
            uqss5=num2cell(ss5l(2:end));
            [S1.Label]=uqss5{:};
            
            %% save segmentation result
            nameFormat1='img_000000%03d_segment_000.jpg';
            fname3= sprintf(nameFormat1,i-1);
            imwrite(ss5, fname3);
            
            %%
            %{
            % test code for label and regionprops is coupled.
            smplist=randi([1,size(S1,1)],1,10);
            disp(smplist);
            for n=1:length(smplist)
                tar=S1(smplist(n)).Label;
                cnt=S1(smplist(n)).Centroid;
                disp(cnt);
                figure; imagesc(ss5==tar);
                hold on;
                scatter(cnt(1), cnt(2), 150,'d');
            end
            %}
            stdInt=zeros(size(S1,1),1); % standard dev of pixel value
            skw=zeros(size(S1,1),1); % skewness of pixel value
            for jj=1:size(S1) % why start from 2
                %for jj=2:size(S1) % why start from 2
                stdInt(jj)=std(double(S1(jj).PixelValues));
                skw(jj)=skewness(double(S1(jj).PixelValues));
            end
            
            %% quantify for resultmat
            S2=regionprops(ss5,B,'MeanIntensity','Centroid','Area',...
                'Perimeter','BoundingBox','PixelIdxList','PixelValues','Eccentricity','EulerNumber');
            delid2=[S2.Area]==0;
            S2(delid2==1)=[];
            
            cent=reshape([S2.Centroid],2,length([S2.Area]));
            T=table(zeros(size([S2.MeanIntensity]'))+i,cent(1,:)',cent(2,:)',...
                [S2.Area]',[S2.Perimeter]',stdInt,[S2.MeanIntensity]', skw);
            T.Properties.VariableNames = {'t','xcoord', 'ycoord', 'area', 'Perimeter',...
                [cName{1},'stdInt'],[cName{1},'meanInt'], 'skewness'};
            segChTable{i}=T;
            
                       
            
            %{
            if i<extraFrameNum
                %exT=[T,ttemp];
                exT=[T,exTab];
                exTable=[exTable;exT];
            end
            %resultTable=[resultTable;T(2:end,:)];
            %}
            %% record other channel's intensity
            intTab_w= cell(1,length(cFormat));

            for ijk=1:length(cFormat)
                fname= sprintf(cFormat{ijk},i-1);
                C=imread(fname);
                C=uint16(C);
                C=double(C);
                S2w=regionprops(ss5,C,'MeanIntensity','PixelIdxList','PixelValues','Area');
                S2w(delid2==1)=[];
                ttemp_w=table([S2w.MeanIntensity]'); % intensity of other channel
                %disp(size(ttemp));
                ttemp_w.Properties.VariableNames={[char(cName(ijk)),'_whole']};
                %sumResultTable{pj,2,nf}=ttemp;
                intTab_w{ijk}=ttemp_w;
            end
            intChTable_w{i}=intTab_w;
            %disp(size(ttemp));
            
            if i<extraFrameNum
                fname= sprintf(cFormat{end},i-1);
                C=imread(fname);
                C=uint16(C);
                C=double(C);
                S11=regionprops(gss2,C,'MeanIntensity','PixelIdxList','PixelValues');
                %S11=regionprops(gather(gss2),C,'MeanIntensity','PixelIdxList','PixelValues');
                S11(delid==1)=[];
                ttemp=table([S11.MeanIntensity]');
                ttemp.Properties.VariableNames=cName(end);
                extChTable{i}=ttemp;
            
            end
            %}
        end
        sumIntTable_w{pj}=intChTable_w;
        sumAreaTable{pj}=areaTable;
        sumSegTable{pj}=segChTable;
        sumExtTable{pj}=extChTable;
        dataPathTable{pj}=foldername;
    end
    
    %% gather quantified data for multiple position into single varialbe
    for pj =1:length(posList{pp}) %position
        
        tmpSegTable=sumSegTable{pj}; % 1st position
        tmpIntTable_w=sumIntTable_w{pj};
        tmpAreaTable=sumAreaTable{pj};
        
        for i=1:length(tmpSegTable) % frame number
            % summlize other channels data to one table
            for k=1:length(tmpIntTable_w{i}) % channel
                if k==1
                    allIntTable= tmpIntTable_w{i}{1};
                else
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
save([datadate,'_segment_',datestr(now,formatOut)],'-v7.3')

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
save([datadate,'_trackinig_',datestr(now,formatOut)],'sumResultTable','sumTrackTable')

%% data cleaning
for i=1:size(sumTrackTable,2)
    tracks=sumTrackTable{i};
    %% cut short trace
    tracks=cutShortTrace(tracks,nf,0.5);
    %% Calc ratiometric value
    tracks=calcRatio(tracks,'FRET_whole','CFP_whole','FC_ratio_whole');
    %% Photo-bleaching correction
    %tracks=pbCorrection(tracks,{'FC_ratio_whole','FC_ratio_whole'});
    sumTrackTable{i}=tracks;
end
save([datadate,'_cleaning_',datestr(now,formatOut)],'sumResultTable','sumTrackTable','dateList','posList')
mvSegmentationResult(dateList, posList,[],[])

%% plot result 
plotParam{1}=sumTrackTable;
plotParam{2}=posList;
plotParam{3}=dateList;
plotParam{4}=genPos;
plotParam{5}='PATH_TO_RESULT_DIRECTORY'; % store directory
plotParam{6}={'FC_ratio'}; %plotType
% left Y-axis
plotParam{7}={'FC_ratio_whole'}; % colum name for left y-axis
plotParam{8}='FRET/CFP Ratio'; % label for left y-axis
% right Y-axis 
plotParam{9}={}; % colum name for right y-axis
plotParam{10}='';% label for right y-axis
% time
plotParam{11}=4; % interval
plotParam{12}='Time (min)'; % label for x-axis
plotParam{13}='all';
plotParam{14}='png';
plotColData(plotParam);

%%
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
plotParam{10}='Standard Deviation of CFP intensity ';% label for right y-axis
% time
plotParam{11}=4; % interval
plotParam{12}='Time (min)'; % label for x-axis
plotParam{13}='all';
plotParam{14}='png';
plotColData(plotParam);
