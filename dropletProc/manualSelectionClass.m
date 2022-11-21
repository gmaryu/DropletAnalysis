classdef manualSelectionClass < handle
    properties (SetAccess = private)
        oriSignal=[];
        procSignal=[];        
    end
    
    properties
        feature;
        yshift;
        manselectClass;
        backClass;
    end
    
    methods
        function obj=manualSelectionClass()
            obj.yshift=0;
            obj.manselectClass=[]; % manualselect
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
            obj.manselectClass=[]; % manualselect
            obj.backClass=[];   
        end
        
        function sig=get.manselectClass(obj)
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
        
        function setRef(obj,ref1,ref2)
             obj.refseq1=ref1;
             obj.refseq2=ref2;
        end
        
        function setOri(obj,input)
            % first peak detection
            obj.oriSignal=input;
            obj.backClass=[];
        end
        
        function centerSig(obj)
            % move signal mean around 0
            obj.yshift=obj.yshift+mean(obj.procSignal);
            obj.procSignal=obj.procSignal-mean(obj.procSignal);
        end
        
        function posiSig(obj)
            % define 
            obj.yshift=obj.yshift+min(obj.procSignal);
            obj.procSignal=obj.procSignal-min(obj.procSignal);
        end
        
        function fillOut(obj)
            % fill NaN value
            obj.procSignal=filloutliers(obj.procSignal,'spline','movmedian',obj.smoothThresh);
        end

        function lowPass(obj)
            % apply low pass filter
            obj.procSignal = filtfilt(obj.lpFilt,obj.procSignal); % Append D zeros to the input data
        end
        
        %{
        function calcBase(obj)
            obj.baseline = filtfilt(obj.baseFilt,obj.procSignal);
        end
        %}
        
        function genEnv(obj)
            % generate an envelope for peak detection. This is critical for
            % peakd etection
            [obj.envu,obj.envl]=envelope(obj.procSignal,obj.peakMaxDist/2,'peak');
            %[obj.envu,obj.envl]=envelope(obj.procSignal, 10, 'rms');
        end
        
        function genPeak(obj)
            % peak detection with processed signal
            [pks,locs,~,p]=findpeaks(obj.procSignal,'MinPeakDistance',obj.peakMinDist);            
            promThresh=max(0,obj.envu-obj.envl)*obj.senThresh+obj.promSens; % prominence threshold array for all timepoints
            pkThresh=promThresh(locs); % prominence threshold array for peak detected timepoints 
            %disp(p);
            %disp(pkThresh);
            % clear peak information if p value (prominence) was smaller than pkthresh
            pks(p<pkThresh)=[];
            locs(p<pkThresh)=[];
            p(p<pkThresh)=[];
            %disp('after')
            %disp(locs);
            %disp(size(locs));
            obj.peakClass=annotationClass(locs,pks+obj.yshift,p);
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom); % property name
        end
        
        function attachPeak(obj)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            % combine processed signal and original signal results
            [pks0,locs0,~,~]=findpeaks(obj.oriSignal);
            
            % assign detection result in genPeak()
            locs=obj.peakClass.time;
            pks=obj.peakClass.value;
            p=obj.peakClass.prom;
            
            
            %{
            for ij=1:length(locs)
                %disp(locs(ij));
                %disp(abs(locs(ij)-locs0));
                %disp(min(abs(locs(ij)-locs0)));
                [v1,id1]=min(abs(locs(ij)-locs0)); % target time point should be selected
                disp(v1);
                msg=['Min: ', num2str(v1), '/id1: ', num2str(id1),'/locs(ij): ',num2str(locs(ij))];
                disp(msg);
                am1=pks0(id1)-pks(ij); % amplitude difference between original peak value and processed peak value
                refamp=median(p(max(1,ij-obj.adaptWindSize):min(length(p),ij+obj.adaptWindSize)));
                % amplicaition range
                %disp(am1/refamp);
                %disp(refamp);
                if am1/refamp<obj.ampThresh(2)&&refamp>obj.ampThresh(1)
                    locs(ij)=locs0(id1);
                    pks(ij)=pks0(id1);
                    %             else
                    %                 locs0(id1)=0;
                end
            end
            % overwrite
            obj.peakClass.time=locs;
            obj.peakClass.value=pks;
            obj.peakClass.prom=p;
            %}
            for ij=1:length(locs)
                %disp(locs(ij));
                %disp(abs(locs(ij)-locs0));
                %disp(min(abs(locs(ij)-locs0)));
                [dif,id1]=min(abs(locs(ij)-locs0)); % target time point should be selected
                if dif > 5
                    disp('this peak wasnt detected in oriSignal');
                    locs(ij)=0;
                else
                    am1=pks0(id1)-pks(ij); % amplitude difference between original peak value and processed peak value
                    refamp=median(p(max(1,ij-obj.adaptWindSize):min(length(p),ij+obj.adaptWindSize)));
                    % amplicaition range
                    %disp(am1/refamp);
                    %disp(refamp);
                    if am1/refamp<obj.ampThresh(2)&&refamp>obj.ampThresh(1)
                        locs(ij)=locs0(id1);
                        pks(ij)=pks0(id1);
                        %             else
                        %                 locs0(id1)=0;
                    end
                end
            end
            locs2=locs(locs~=0);
            pks2=pks(locs~=0);
            p2=p(locs~=0);
            % overwrite
            obj.peakClass.time=locs2;
            obj.peakClass.value=pks2;
            obj.peakClass.prom=p2;
        end
        
        function genTrough(obj)
            % detection of trough (bottom peak)
            locP=obj.peakClass.time % location of peak
            loc=zeros(length(locP)-1,1);
            value=zeros(length(locP)-1,1);
            prom=zeros(length(locP)-1,1);
            counter=zeros(size(obj.oriSignal));
            counter(locP)=1;
            counter=cumsum(counter); % value is the serial number of peaks
            for ii=1:length(locP)-1
                [~,idx]=min(obj.oriSignal(counter==ii)); % min index of specific counter value 
                loc(ii)=idx+find(counter==ii, 1)-1; % convert relative timepoint to absolute timepoint
                value(ii)=obj.oriSignal(loc(ii)); % get real value of trough point
                % this algorithm is not robust for some data that doesn't
                % have clear bottom value
            end
            obj.troughClass=annotationClass(loc,value,prom);
            
        end
        
        function showSignal(obj)
            % plot original signal, processed signal for peak detection,
            % envelope upper and lower score.
            figure(1)
            plot(obj.oriSignal)
            hold on
            plot(obj.procSignal+obj.yshift)
            plot(obj.envu+obj.yshift)
            plot(obj.envl+obj.yshift)
            hold off
        end
        
        function showPeak(obj)
            % show peaks only
            figure(1);clf;hold on
            plot(obj.oriSignal,'r')
            % plot(obj.procSignal+obj.yshift,'.-c')
            scatter(obj.peakClass.time,obj.peakClass.value,60,'bv')
            plot(obj.envu+obj.yshift)
            plot(obj.envl+obj.yshift)
            % text(obj.peakClass.time,obj.peakClass.value+2,...
            %       cellstr(num2str((1:length(obj.peakClass.value))')))
            hold off
        end
        
        function userTPSelect(obj,t)
            tt=round(t);
            vv=obj.oriSignal(tt);
            obj.manselectClass=annotationClass2(t,vv);
            figure(1)
            hold on
            scatter(tt,vv,60,'bv')
            hold off
            obj.manselectClass.addAnnoTime(tt,vv);
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
        
        function checkPeak(obj)
            obj.backClass=annotationClass(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom);
            tt=obj.peakClass.time;
            vv=obj.peakClass.value;
            pp=obj.peakClass.prom;
            pm=pp; % peak max?
            for kk=2:length(pp)
                pm(kk)=max(pp(1:kk)); % update maximum peak value?
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
        
        function [drop,osci,raw]=output(obj,n)
            % generate a table for summilize peak and trough data
            peakId=zeros(size(obj.oriSignal));
            troughId=zeros(size(obj.oriSignal));
            peakId(obj.peakClass.time)=1;
            troughId(obj.troughClass.time)=1;
            peakCount=cumsum(peakId);
            troughCount=cumsum(troughId);
            rId=1:length(peakId);
            raw=table(rId',obj.oriSignal,obj.procSignal,obj.envu,obj.envl,...
                        peakId,troughId,peakCount,troughCount,zeros(size(rId'))+n);
            raw=[raw,obj.feature];
            dId=1:length(obj.peakClass.time)
            ampLeft=[NaN;obj.peakClass.value(2:end)-obj.troughClass.value]
            ampRight=[NaN;obj.peakClass.value(1:end-1)-obj.troughClass.value]
            timeRight=[NaN;-obj.peakClass.time(1:end-1)+obj.troughClass.time]
            timeLeft=[NaN;obj.peakClass.time(2:end)-obj.troughClass.time]         
            periodPeak=[NaN;obj.peakClass.time(2:end)-obj.peakClass.time(1:end-1)]
            periodTrough=[NaN;obj.troughClass.time(2:end)-obj.troughClass.time(1:end-1);NaN]
            try
%                 disp(obj.peakClass.time);
%                 disp(obj.peakClass.value);
%                 disp(obj.peakClass.prom);
%                 disp([nan;obj.troughClass.time]);
%                 disp([nan;obj.troughClass.value]);
%                 disp([nan;obj.troughClass.prom]);
%                 disp(ampLeft);
%                 disp(ampRight);
%                 disp(timeLeft);
%                 disp(timeRight);
%                 disp(periodPeak);
%                 disp(periodTrough);
%                 disp(dId');
%                 disp(zeros(size(dId'))+n);
                osci=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,...
                        [nan;obj.troughClass.time],[nan;obj.troughClass.value],[nan;obj.troughClass.prom],...
                        ampLeft,ampRight,timeLeft,timeRight,periodPeak,periodTrough,dId',zeros(size(dId'))+n);
            catch
%                 disp(obj.peakClass.time);
%                 disp(obj.peakClass.value);
%                 disp(obj.peakClass.prom);
%                 disp([nan;obj.troughClass.time]);
%                 disp([nan;obj.troughClass.value]);
%                 disp([nan;obj.troughClass.prom]);
%                 disp(ampLeft);
%                 disp(ampRight);
%                 disp(timeLeft);
%                 disp(timeRight);
%                 disp(periodPeak);
%                 disp(size(periodTrough));
%                 disp(dId');
%                 disp(zeros(size(dId'))+n);
                osci=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,...
                            [nan;obj.troughClass.time],[nan;obj.troughClass.value],[nan;obj.troughClass.prom],...
                            ampLeft,ampRight,timeLeft,timeRight,periodPeak,periodTrough(1:end-1),dId',zeros(size(dId'))+n);
                obj
            end

            raw.Properties.VariableNames={'frameID','oriSignal','procSignal','envu','envl',...
                'peakId','troughId','peakCount','troughCount','dropID','frame','madSig',...
                'areaSig','xcoorSig','ycoorSig'};
            osci.Properties.VariableNames={'peakTime','peakValue','peakProm','troughTime',...
                'troughValue','troughProm','ampLeft','ampRight','timeLeft','timeRight',...
                'periodPeak','periodTrough','cycleID','dropID'};
            
            drop=table(nanmean(periodPeak),nanmean(periodTrough),mad(periodPeak),mad(periodTrough),...
                nanmean(ampLeft),nanmean(ampRight),mad(ampLeft),mad(ampRight),n);
            drop.Properties.VariableNames={'meanPerPeak','meanPerTrough','madPerPeak','madPerTrough',...
                'meanAmpLeft','meanAmpRight','madAmpLeft','madAmpRight','dropID'};
        end
        
    end
end