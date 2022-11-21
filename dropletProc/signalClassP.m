classdef signalClassP < handle
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
        peakClass;
        backClass;
    end
    
    methods
        function obj=signalClassP(sen,prom,smooth,peak,...
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
            obj.oriPeakClass=[];
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
        end
        
        function clearObj(obj)
            obj.yshift=0;
            obj.peakClass=[];
            obj.oriSignal=[];
            obj.procSignal=[];
            obj.oriPeakClass=[];
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
            % first peak detection
            obj.oriSignal=input;
            obj.procSignal=input;
            [pks,locs,~,p]=findpeaks(obj.procSignal);
            obj.oriPeakClass=annotationClass(locs,pks,p);
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
            %[obj.envu,obj.envl]=envelope(obj.procSignal);
        end
        
        function genPeak(obj)
            % peak detection with processed signal
            [pks,locs,~,p]=findpeaks(obj.procSignal,'MinPeakDistance',obj.peakMinDist);
            promThresh=max(0,obj.envu-obj.envl)*obj.senThresh+obj.promSens; % prominence threshold array for all timepoints
            pkThresh=promThresh(locs); % prominence threshold array for peak detected timepoints 
            % clear peak information if p value (prominence) was smaller than pkthresh
            pks(p<pkThresh)=[];
            locs(p<pkThresh)=[];
            p(p<pkThresh)=[];
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
            
            for ij=1:length(locs)
                [~,id1]=min(abs(locs(ij)-locs0)); % target time point should be selected
                am1=pks0(id1)-pks(ij);
                refamp=median(p(max(1,ij-obj.adaptWindSize):min(length(p),ij+obj.adaptWindSize)));
                % amplicaition range
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
        end
        
        
        function showSignal(obj)
            % plot original signal, processed signal for peak detection,
            % envelope upper and lower score.
            figure(10)
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

        
        function [drop,osci,raw]=output(obj,n)
            disp(n)
            % generate a table for summilize peak and trough data
            peakId=zeros(size(obj.oriSignal));
            peakId(obj.peakClass.time)=1;
            peakCount=cumsum(peakId);
            rId=1:length(peakId);
            raw=table(rId',obj.oriSignal,obj.procSignal,obj.envu,obj.envl,...
                        peakId,peakCount,zeros(size(rId'))+n);
            raw=[raw,obj.feature];
            dId=1:length(obj.peakClass.time);
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
                osci=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,dId',zeros(size(dId'))+n);
                disp('This works!')
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
                osci=table(obj.peakClass.time,obj.peakClass.value,obj.peakClass.prom,dId',zeros(size(dId'))+n);
            end

            raw.Properties.VariableNames={'frameID','oriSignal','procSignal','envu','envl',...
                'peakId','peakCount','dropID','frame','madSig',...
                'areaSig','xcoorSig','ycoorSig'};
            osci.Properties.VariableNames={'peakTime','peakValue','peakProm','cycleID','dropID'};
            
            drop=table(n);
            drop.Properties.VariableNames={'dropID'};
        end
        
    end
end