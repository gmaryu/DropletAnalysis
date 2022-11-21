function [tracks,trackArray] = CellTrackerNostat(feat,refCh,parset)
% cell tracker in Live mode
tic

if ~exist ('parset','var')
    par=zeros(13,1);
    par(1)=30;
    par(2)=10;
    par(3)=std(log10(feat.area)); %log scale
    par(4)=std(feat.mCherry);
    par(5:8)=1;
    par(9:12)=0;
else
    par=parset;
end

maxDisp=par(1);
%Parse inputs
t=feat.t;
xbound=max(feat.xcoord);
ybound=max(feat.ycoord);
[tt,begins]=unique(t); %tt: start of time point
ends=circshift(begins,[-1 0])-1;
ends(length(ends))=length(t); %ends: ends of time point
Nf=numel(tt); % number of time point
trackArray=zeros(size(feat,1),1); % track lable
segState=zeros(size(feat,1),1); % segmentation overlap label

if Nf < 2
    error(['Sorry, found too few files'])
end

%First frame
ind=begins(1):ends(1);
nparticles = numel(ind);
tracks = repmat(struct('Len',[],'Leaf',[],'Stem',[],'Start',[],'Feat',[],'State',[],'Statis',[]), ...
    nparticles,1);
for ii = 1:nparticles
    tracks(ii) = struct('Len',1,'Leaf',0,'Stem',0,'Start',1,'Feat',feat(ind(ii),:),'State',0,'Statis',[]);
end
trackArray(1:nparticles)=1:nparticles;
KDTold=KDTreeSearcher([feat.xcoord(ind),feat.ycoord(ind)]);

%Loop over frames
for t = 2:Nf
%     if t==57
%         aaaa=1;
%     end
    ind=begins(t):ends(t); % new frame
    ind0=begins(t-1):ends(t-1); % this frame
    nfr1 = numel(ind); % number of segments in new frame
    nfr0 = numel(ind0);% number of segments in this frame
    if nfr1==0
        warning('MATLAB:CellTracker:noCells', ...
            ['Found no cells in frame ' num2str(t) '.']);
    end
    %feature of new frame
    %nextF=feat{ind,1:end-1};
    %nextF=[feat.t(ind),feat.xcoord(ind),feat.ycoord(ind),feat.area(ind),feat.FRET(ind),feat.CFP(ind)];
    %nextF=[feat.t(ind),feat.xcoord(ind),feat.ycoord(ind),feat.area(ind),feat.mCherrymeanInt(ind),feat.CFP(ind)];
    nextF= [feat.t(ind),feat.xcoord(ind),feat.ycoord(ind),feat.area(ind),eval(['feat.',refCh,'meanInt(ind)'])];
    %this frame
    %thisF = [feat.t(ind0),feat.xcoord(ind0),feat.ycoord(ind0),feat.area(ind0),feat.FRET(ind0),feat.CFP(ind0)];
    %thisF = [feat.t(ind0),feat.xcoord(ind0),feat.ycoord(ind0),feat.area(ind0),feat.mCherrymeanInt(ind0),feat.CFP(ind0)];
    thisF=[feat.t(ind0),feat.xcoord(ind0),feat.ycoord(ind0),feat.area(ind0),eval(['feat.',refCh,'meanInt(ind0)'])];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% frame before
%     prior = zeros(size(thisF));
% calculate spped
%     for ii = 1:size(thisF,1)
%         tr = tracks(trackArray(ind0(ii)));
%         if tr.Len > 1
%             prior(ii,:) = tr.Feat(end-1,:);
%         else
%             prior(ii,:) = thisF(ii,:);
%         end
%     end
% 
%     % estimate a velocity for each particle in fr0
%     velocity = thisF(:,[2,3]) - prior(:,[2,3]);
%     % and use kinematics to estimate a future position
%     estimate = thisF(:,[2,3]) + velocity;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % define cost and link arrays
    KDTnew=KDTreeSearcher(nextF(:,2:3));
    if nfr1>0
        % loop over active tracks
        moveIdx=[];
        boundIdx=[];
        for ii = 1:size(thisF,1)           
            % compare the new and old position
            devFlag=[];
            idxNextF=rangesearch(KDTnew,thisF(ii,2:3),maxDisp);
            idxNextF=cell2mat(idxNextF);
            moveIdx=[moveIdx;[zeros(length(idxNextF),1)+ii,idxNextF']];
            %at the boundary?
            boundDist=min([thisF(ii,2),thisF(ii,3),abs(xbound-thisF(ii,2)),abs(ybound-thisF(ii,3))]);
            boundFlag=boundDist<maxDisp;
            boundIdx=[boundIdx;boundFlag*boundDist];

            % going division?
        end
         %relation mat
        Cmat=[];
        if ~isempty(moveIdx)
        Cmat=[Cmat;sparse(repmat(1:size(moveIdx,1),[1,2]),...
            [moveIdx(:,1);nfr0+moveIdx(:,2)],1,size(moveIdx,1),nfr1+nfr0)];
        end
        Cmat=[Cmat;sparse(1:nfr0,1:nfr0,1,nfr0,nfr1+nfr0)];
        %prob mat       
        Pmat=zeros(size(Cmat,1),1);
        count=0;
        %calc p
        % normal
        for jj=1:size(moveIdx,1)
            count=count+1;
            Pmat(count)=calcFeatProb(thisF(moveIdx(jj,1),:),nextF(moveIdx(jj,2),:),'move',...
                par(2:end)',1);
        end
        for jj=1:nfr0
            if boundIdx(jj)>0
                count =count+1;
                Pmat(count)=calcFeatProb(boundIdx(jj),[],'moveout',...
                    par(2:end)',1);
            else
                count=count+1;
                Pmat(count)=calcFeatProb([],[],'disappear',...
                    par(2:end)',1);
            end
        end
        %
        boundIdx=zeros(size(nextF,1),1);
        for ii=1:size(nextF,1)
            boundDist=min([nextF(ii,2),nextF(ii,3),abs(xbound-nextF(ii,2)),abs(ybound-nextF(ii,3))]);
            boundFlag=boundDist<maxDisp;
            boundIdx(ii)=boundFlag*boundDist;           
        end
        Cmat2=[];

        Cmat2=[Cmat2;sparse(1:nfr1,nfr0+(1:nfr1),1,nfr1,nfr1+nfr0)];
        Pmat2=zeros(size(Cmat2,1),1);
        count=0;
        for jj=1:nfr1
            if boundIdx(jj)>0
                count =count+1;
                Pmat2(count)=calcFeatProb(boundIdx(jj),[],'movein',...
                    par(2:end)',1);
            else
                count=count+1;
                Pmat2(count)=calcFeatProb([],[],'appear',...
                    par(2:end)',1);
            end
        end
        Pmat=[Pmat;Pmat2];
        Cmat=[Cmat;Cmat2];
%         if adaptMaxDisp&&maxDisp>maxDispMin
%             if(length(Pmat)>2000)
%                 maxDisp=maxDisp*0.9;
%                 disp(['Mas Disp Change to', num2str(maxDisp)])
%             end
%         end

        options = optimoptions(@intlinprog,'Display','off');
        link=intlinprog(-(Pmat),1:length(Pmat),[],[],Cmat(:,1:end)',...
            ones(size(Cmat,2),1),zeros(length(Pmat),1),ones(length(Pmat),1),options);
        link=Cmat(abs(link-1)<0.1,:);
        
        for ii=1:size(link,1)
            src=find(link(ii,[1:size(thisF,1)])==1);
            N_src=length(src);
            dest=find(link(ii,[size(thisF,1)+1:size(thisF,1)+size(nextF,1)])==1);
            N_dest=length(dest);
            if N_src==1
                if N_dest==1
                        trackArray(ind(dest))=trackArray(ind0(src));
                        segState(ind(dest))=segState(ind0(src));
                        tracks(trackArray(ind0(src))).Feat(end+1,:)=...
                            feat(ind(dest(N_dest)),:);
                        tracks(trackArray(ind0(src))).Len=tracks(trackArray(ind0(src))).Len+1;     
                end
                if N_dest>1 
                    warning('MATLAB:CellTracker:linkError',...
                    ['More than one dest found!']);
                end
            elseif N_src==0
                if N_dest~=1
                    warning('MATLAB:CellTracker:linkError', ...
                        ['Cannot find link!']);
                else        
                    trackArray(ind(dest(1)))=length(tracks)+1;
                    tracks(end+1)=struct('Len',1,'Leaf',0,'Stem',0,'Start',t,'Feat',feat(ind(dest(1)),:),'State',0,'Statis',[]); 
                end
            end            
        end        
    end
    for i=1:length(tracks)
        if tracks(i).State(1)
            if tracks(i).State(3)==1&&tracks(i).State(1)~=0
%                if tracks(tracks(i).State(1)).Feat{end,1}==t
                    tracks(i).Feat(end+1,:)=tracks(tracks(i).State(1)).Feat(end,:);
                    tracks(i).Len= tracks(i).Len+1;
%                else
%                    tracks(i).State(3)=0;
%                end
            end
        end
    end
    disp([num2str(t),'processed'])
    KDTold=KDTnew;
end      
toc
end

