classdef signalClass < handle
    properties (SetAccess = private)
        oriSignal=[];
        procSignal=[];
        baseline=[];
        refseq1=[];
        refseq2=[];
        envu=[];
        envl=[];
        
    end
    properties
        feature;
        yshift;
        senThresh;
        promSens;
        smoothThresh;
        peakMinDist;
        peakMaxDist;
        filtThresh;
        lpFilt;
        adaptWindSize;
        ampThresh;
        % baseFilt;
        maxDelete=5;
        scale=3;
        % adValue=10;
        oriPeakClass;
        oriTroughClass;
        peakClass;
        troughClass;
        backClass;
    end
    
    methods
        function obj=signalClass(sen,prom,smooth,peak,...
                filt,lp,peakM,adws,athresh)
            obj.senThresh=sen;
            obj.promSens=prom;
            obj.smoothThresh=smooth;
            obj.peakMinDist=peak;
            obj.filtThresh=filt;
            obj.lpFilt=lp;
            % obj.baseFilt=base;
            obj.peakMaxDist=peakM;
            obj.adaptWindSize=adws;
            obj.ampThresh=athresh;
            obj.yshift=0;
            obj.peakClass=[];
            obj.troughClass=[];
            obj.oriPeakClass=[];
            obj.oriTroughClass=[];
            obj.backClass=[];
            
        end
        function sigPos=updateObj(obj,raw,osci,dropId,posIdx)
            sigPos=(raw.dropID==dropId)&(raw.posIdx==posIdx);
            tempRaw=raw(sigPos,:);
            tempOsci=osci(sigPos,:);
            obj.oriSignal=tempRaw.oriSignal;
            obj.procSignal=tempRaw.procSignal;
            obj.yshift=mean(obj.oriSignal)-mean(obj.procSignal);
            obj.envu=tempRaw.envu;
            obj.envl=tempRaw.envl;
            obj.peakClass=annotationClass(tempOsci.peakTime,tempOsci.peakValue,tempOsci.peakProm);
            obj.troughClass=annotationClass(tempOsci.troughTime,tempOsci.troughValue,tempOsci.troughProm);
            
        end
        function clearObj(obj)
            obj.yshift=0;
            obj.peakClass=[];
            obj.troughClass=[];
            obj.oriSignal=[];
            obj.procSignal=[];
            obj.oriPeakClass=[];
            obj.oriTroughClass=[];
            obj.baseline=[];
            obj.envu=[];
            obj.envl=[];
            obj.backClass=[];
        end
        function sig=get.oriSignal(obj)
            sig=obj.oriSignal;
        end
        function sig=get.procSignal(obj)
            sig=obj.procSignal;
        end
        function sig=get.feature(obj)
            sig=obj.feature;
        end
        
        function sig=get.yshift(obj)
            sig=obj.yshift;
        end
        function sig=get.refseq1(obj)
            sig=obj.refseq1;
        end
        function sig=get.refseq2(obj)
            sig=obj.refseq2;
        end
        function setRef(obj,ref1,ref2)
             obj.refseq1=ref1;
             obj.refseq2=ref2;
        end
        function setOri(obj,input)
            obj.oriSignal=input;
            %disp(max(obj.oriSignal));
            obj.procSignal=input;
            %disp(max(obj.procSignal));
            [pks,locs,~,p]=findpeaks(obj.procSignal);
            obj.oriPeakClass=annotationClass(locs,pks,p);
            [pks,locs,~,p]=findpeaks(-obj.procSignal);
            obj.oriTroughClass=annotationClass(locs,-pks,p);
            obj.backClass=[];
            % so many peaks generated in this section
        end
        
        function centerSig(obj)
            obj.yshift=obj.yshift+mean(obj.procSignal);
            obj.procSignal=obj.procSignal-mean(obj.procSignal);
        end
        function posiSig(obj)
            obj.yshift=obj.yshift+min(obj.procSignal);
            obj.procSignal=obj.procSignal-min(obj.procSignal);
        end
        function fillOut(obj)
            obj.procSignal=filloutliers(obj.procSignal,'spline','movmedian',obj.smoothThresh);
        end
        %ym=smooth(ym,12,'sgolay');
        function lowPass(obj)
            obj.procSignal = filtfilt(obj.lpFilt,obj.procSignal); % Append D zeros to the input data
        end
        %         function calcBase(obj)
        %             obj.baseline = filtfilt(obj.baseFilt,obj.procSignal);
        %         end
        function genEnv(obj)
            %[obj.envu,obj.envl]=envelope(obj.procSignal,obj.peakMaxDist/2,'peak');
            [obj.envu,obj.envl]=envelope(obj.procSignal);
            
        end
        function genPeak(obj)
            [pks,locs,~,p]=findpeaks(obj.procSignal,'MinPeakDistance',obj.peakMinDist);
            %disp(size(pks));
            promThresh=max(0,obj.envu-obj.envl)*obj.senThresh+obj.promSens; % prominence threshold array for all timepoints
            %disp(promThresh);
            pkThresh=promThresh(locs); % prominence threshold array for peak detected timepoints 
            %disp(pkThresh);
            pks(p<pkThresh)=[];locs(p<pkThresh)=[];p(p<pkThresh)=[]; % if p value (prominence) was smaller than pkthresh, delete its data
            obj.peakClass=annotationClass(locs,pks+obj.yshift,p);
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom); % property name
        end
        function attachPeak(obj)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            % combine processed signal and original signal results?
            [pks0,locs0,~,~]=findpeaks(obj.oriSignal);
            
            % assign detection result in genPeak()
            locs=obj.peakClass.time;
            pks=obj.peakClass.value;
            p=obj.peakClass.prom;
            
            for ij=1:length(locs)
                [~,id1]=min(abs(locs(ij)-locs0)); % target time point should be selected
                am1=pks0(id1)-pks(ij);
                refamp=median(p(max(1,ij-obj.adaptWindSize):min(length(p),ij+obj.adaptWindSize)));
                if am1/refamp<obj.ampThresh(2)&&refamp>obj.ampThresh(1)
                    locs(ij)=locs0(id1);
                    pks(ij)=pks0(id1);
                    %             else
                    %                 locs0(id1)=0;
                end
            end
            obj.peakClass.time=locs;
            obj.peakClass.value=pks;
            obj.peakClass.prom=p;
        end
        function genTrough(obj)
            locP=obj.peakClass.time;
            loc=zeros(length(locP)-1,1);
            value=zeros(length(locP)-1,1);
            prom=zeros(length(locP)-1,1);
            counter=zeros(size(obj.oriSignal));
            counter(locP)=1;
            counter=cumsum(counter);
            for ii=1:length(locP)-1
                [~,idx]=min(obj.oriSignal(counter==ii));
                loc(ii)=idx+find(counter==ii, 1 )-1;
                value(ii)=obj.oriSignal(loc(ii));
            end
            obj.troughClass=annotationClass(loc,value,prom);
            
        end
        function showSignal(obj)
            figure(10)
            plot(obj.oriSignal)
            hold on
            plot(obj.procSignal+obj.yshift)
            plot(obj.envu+obj.yshift)
            plot(obj.envl+obj.yshift)
            hold off
        end
        function showPeak(obj)
            figure(1)
            clf
%             plot(obj.procSignal+obj.yshift,'.-c')
            hold on
            plot(obj.oriSignal,'r')
            scatter(obj.peakClass.time,obj.peakClass.value,60,'bv')
            %             text(obj.peakClass.time,obj.peakClass.value+2,...
            %                 cellstr(num2str((1:length(obj.peakClass.value))')))
            hold off
        end
        function showTrough(obj)
            figure(1)
            clf
           % plot(obj.procSignal+obj.yshift,'.-c')
            hold on
            plot(obj.oriSignal,'r')
            scatter(obj.peakClass.time,obj.peakClass.value,60,'cv')
            scatter(obj.troughClass.time,obj.troughClass.value,60,'b^')
            %             text(obj.troughClass.time,obj.troughClass.value+2,...
            %                 cellstr(num2str((1:length(obj.troughClass.value))')))
            hold off
        end
        function userAddPeak(obj,t)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            [~,idx]=min(abs(obj.oriPeakClass.time-t));
            tt=obj.oriPeakClass.time(idx);
            if ~isempty(find(obj.peakClass.time==tt, 1))
                warning('peak already in list')
                return;
            end
            vv=obj.oriPeakClass.value(idx);
            pp=obj.oriPeakClass.prom(idx);
            figure(1)
            hold on
            scatter(tt,vv,60,'bv')
            hold off
            obj.peakClass.addAnnoTime(tt,vv,pp);
        end
        function userAddTrough(obj,t)
            obj.backClass=annotationClass(obj.troughClass.time,obj.troughClass.value,obj.troughClass.prom);
            [~,idx]=min(abs(obj.oriTroughClass.time-t));
            tt=obj.oriTroughClass.time(idx);
            if ~isempty(find(obj.troughClass.time==tt, 1))
                warning('trough already in list')
                return;
            end
            vv=obj.oriTroughClass.value(idx);
            pp=obj.oriTroughClass.prom(idx);
            figure(1)
            hold on
            scatter(tt,vv,60,'b^')
            hold off
            obj.troughClass.addAnnoTime(tt,vv,pp);
        end
        function userDelPeak(obj,t)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            [~,idx]=min(abs(obj.peakClass.time-t));
            tt=obj.peakClass.time(idx);
            vv=obj.peakClass.value(idx);
            figure(1)
            hold on
            scatter(tt,vv,60,'wv')
            hold off
            obj.peakClass.delAnnoTime(tt);
        end
        function userDelTrough(obj,t)
            obj.backClass=annotationClass(obj.troughClass.time,obj.troughClass.value,obj.troughClass.prom);
            [~,idx]=min(abs(obj.peakClass.time-t));
            tt=obj.troughClass.time(idx);
            vv=obj.troughClass.value(idx);
            figure(1)
            hold on
            scatter(tt,vv,60,'w^')
            hold off
            obj.troughClass.delAnnoTime(tt);
        end
        
        function checkPeak(obj)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            tt=obj.peakClass.time;
            vv=obj.peakClass.value;
            pp=obj.peakClass.prom;
            pm=pp;
            for kk=2:length(pp)
                pm(kk)=max(pp(1:kk));
            end
            vt=obj.troughClass.value;
            ttrough=obj.troughClass.time;
            pl=vv-[min(obj.oriSignal(1:ttrough(1)));vt];
            pr=vv-[vt;min(obj.oriSignal(ttrough(end):end))];
            
            ppp=(pl+pr)/prctile((pl+pr),75);
            ttn=tt/prctile(tt,75);
            valMat=[ppp,ttn];
            stats.mad_s=inf;
            for ik=1:10
                try
                    [~,tstat]=robustfit(valMat(:,1),valMat(:,2),'bisquare',2.^(ik-7));
                catch
                    obj
                end
                if tstat.mad_s <stats.mad_s
                    stats=tstat;
                end
            end
            %
            %             figure(2)
            %             scatter(valMat(:,1),valMat(:,2))
            %             hold on
            %             plot(0:0.1:8,b(1)+b(2)*(0:0.1:8))
            %             hold off
            %
            %             figure(3)
            %             plot(abs(log10(pl./pr)))
            %             hold off
            
            balanceIdx=abs(log10(pl./pr));
            delIdx=(abs(stats.resid )>obj.filtThresh(1));
            delIdx=delIdx|(balanceIdx>obj.filtThresh(2));
            delIdx=delIdx|(pm./pp>obj.filtThresh(3));
            
            if ~isempty(find(delIdx(1:3)==1, 1))
                delIdx(1:3)=1;
            end
            pp(delIdx)=[];tt(delIdx)=[];vv(delIdx)=[];
            obj.peakClass.time=tt;
            obj.peakClass.value=vv;
            obj.peakClass.prom=pp;
            obj.genTrough();
        end
        
        function [osci,drop,raw]=output(obj,n)
            peakId=zeros(size(obj.oriSignal));
            %troughId=zeros(size(obj.oriSignal));
            peakId(obj.peakClass.time)=1;
            %troughId(obj.troughClass.time)=1;
            peakCount=cumsum(peakId);
            %troughCount=cumsum(troughId);
            rId=1:length(peakId);
            
            raw=table(rId',obj.oriSignal,obj.procSignal,obj.envu,obj.envl,...
                peakId,peakCount,zeros(size(rId'))+n);
            %raw=table(rId',obj.oriSignal,obj.procSignal,obj.envu,obj.envl,...
            %    peakId,troughId,peakCount,troughCount,zeros(size(rId'))+n);
            raw=[raw,obj.feature];
            dId=1:length(obj.peakClass.time);
            ampLeft=[NaN;obj.peakClass.value(2:end)];
            ampRight=[NaN;obj.peakClass.value(1:end-1)];
            timeRight=[NaN;-obj.peakClass.time(1:end-1)];
            timeLeft=[NaN;obj.peakClass.time(2:end)];
            %{
            ampLeft=[NaN;obj.peakClass.value(2:end)-obj.troughClass.value];
            ampRight=[NaN;obj.peakClass.value(1:end-1)-obj.troughClass.value];
            timeRight=[NaN;-obj.peakClass.time(1:end-1)+obj.troughClass.time];
            timeLeft=[NaN;obj.peakClass.time(2:end)-obj.troughClass.time];
            %}
            
            periodPeak=[NaN;obj.peakClass.time(2:end)-obj.peakClass.time(1:end-1)];
            %periodTrough=[NaN;obj.troughClass.time(2:end)-obj.troughClass.time(1:end-1);NaN];
            try
                %{
                drop=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,...
                    ampLeft,ampRight,timeLeft,timeRight,periodPeak,dId',zeros(size(dId'))+n);
                %}
                drop=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,...
                    [nan;obj.troughClass.time],[nan;obj.troughClass.value],[nan;obj.troughClass.prom],...
                    ampLeft,ampRight,timeLeft,timeRight,periodPeak,periodTrough,dId',zeros(size(dId'))+n);
                    %}
            catch
                obj
            end
            raw.Properties.VariableNames={'frameID','oriSignal','procSignal','envu','envl',...
                'peakId','peakCount','dropID','frame','madSig',...
                'areaSig','xcoorSig','ycoorSig'};
            drop.Properties.VariableNames={'peakTime','peakValue','peakProm',...
                'troughValue','ampLeft','ampRight','timeLeft','timeRight'...
                ,'periodPeak','cycleID','dropID'};
            
            osci=table(nanmean(periodPeak),mad(periodPeak),...
                nanmean(ampLeft),nanmean(ampRight),mad(ampLeft),mad(ampRight),n);
            osci.Properties.VariableNames={'meanPerPeak','madPerPeak',...
                'meanAmpLeft','meanAmpRight','madAmpLeft','madAmpRight','dropID'};
            %{
                            raw.Properties.VariableNames={'frameID','oriSignal','procSignal','envu','envl',...
                'peakId','troughId','peakCount','troughCount','dropID','frame','madSig',...
                'areaSig','xcoorSig','ycoorSig'};
            drop.Properties.VariableNames={'peakTime','peakValue','peakProm','troughTime',...
                'troughValue','troughProm','ampLeft','ampRight','timeLeft','timeRight'...
                ,'periodPeak','periodTrough','cycleID','dropID'};
            
            osci=table(nanmean(periodPeak),nanmean(periodTrough),mad(periodPeak),mad(periodTrough),...
                nanmean(ampLeft),nanmean(ampRight),mad(ampLeft),mad(ampRight),n);
            osci.Properties.VariableNames={'meanPerPeak','meanPerTrough','madPerPeak','madPerTrough',...
                'meanAmpLeft','meanAmpRight','madAmpLeft','madAmpRight','dropID'};
                %}
        end
        
    end
end