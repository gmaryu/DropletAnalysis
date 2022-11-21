classdef annotationClass2<handle
    properties
        time;
        value;
    end
    methods
        function obj=annotationClass2(t,v)
            obj.time=t;
            obj.value=v;
        end
        function delAnnoIndex(obj,in)%delete base on index
            obj.time(in)=[];
            obj.value(in)=[];
        end
        function delAnnoTime(obj,t)%delete all points close to t vector
            in=zeros(size(t),'logical');
            for ii=1:length(t)
                [~,tIdx]=min(abs(t(ii)-obj.time));
                in(tIdx)=1;
            end
            obj.time(in)=[];
            obj.value(in)=[];
        end
        function addAnnoTime(obj,t,v)%add one point
            indx=length(find(t>obj.time))+1;
            obj.time=[obj.time(1:indx-1);t;obj.time(indx:end)];
            obj.value=[obj.value(1:indx-1);v;obj.value(indx:end)];
        end
        function [to,vo]=getAnnoTime(obj,t)
            in=zeros(size(t),'logical');
            
            for ii=1:length(t)
                [~,tIdx]=min(abs(t(ii)-obj.time));
                in(tIdx)=1;
            end
            to=obj.time(in);
            vo=obj.value(in);
        end
    end
end
