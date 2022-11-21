classdef annotationClass<handle
    properties
        time;
        value;
        prom;
    end
    methods
        function obj=annotationClass(t,v,p)
            obj.time=t;
            obj.value=v;
            obj.prom=p;
        end
        function delAnnoIndex(obj,in)%delete base on index
            obj.time(in)=[];
            obj.value(in)=[];
            obj.prom(in)=[];
        end
        function delAnnoTime(obj,t)%delete all points close to t vector
            in=zeros(size(t),'logical');
            for ii=1:length(t)
                [~,tIdx]=min(abs(t(ii)-obj.time));
                in(tIdx)=1;
            end
            obj.time(in)=[];
            obj.value(in)=[];
            obj.prom(in)=[];
        end
        function addAnnoTime(obj,t,v,p)%add one point
            indx=length(find(t>obj.time))+1;
            obj.time=[obj.time(1:indx-1);t;obj.time(indx:end)];
            obj.value=[obj.value(1:indx-1);v;obj.value(indx:end)];
            obj.prom=[obj.prom(1:indx-1);p;obj.prom(indx:end)];
        end
        function [to,vo,po]=getAnnoTime(obj,t)
            in=zeros(size(t),'logical');
            
            for ii=1:length(t)
                [~,tIdx]=min(abs(t(ii)-obj.time));
                in(tIdx)=1;
            end
            to=obj.time(in);
            vo=obj.value(in);
            po=obj.prom(in);
        end
    end
end
