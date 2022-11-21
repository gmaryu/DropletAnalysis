function o=estimateParameter(seg,nChannel)
resmat=[];
mpdm=[];
for i=2:max(seg{:,1})
    idx=seg{:,1}==i;
    idx0=seg{:,1}==i-1;
    
    mat=seg{idx,2:3};
    mat0=seg{idx0,2:3};
    distmat=pdist2(mat,mat0);
    num=round(sqrt(length(mat)*length(mat0)));
    distmat=min(distmat);
 %   distmat=sort(distmat(:),'ascend');
 %   mpdm(i)=mean(distmat);
    mpdm(i)=quantile(distmat,0.9);
    
    vpdm(i)=std(distmat);
%    resmat=[resmat;[distmat(1:num),zeros(num,1)+i]];
end
% range=[-max(abs(seg{:,2}))+100,max(abs(seg{:,2}))+100];
% [p1,x1]=ksdensity(seg{:,2},'support',[0 range(2)+1e-6]);
% range=[-max(abs(x1)),max(abs(x1))];
% ninterp=200;
% px1=interp1([range(1)-1e-6,x1,range(2)+1e-6],[0,p1,0],linspace(range(1),range(2),ninterp));
% px2=interp1([range(1)-1e-6,-x1,range(2)+1e-6],[0,p1(end:-1:1),0],linspace(range(1),range(2),ninterp));
% px1x2=conv(px1,px2);
% range=[-max(abs(seg{:,2}))+100,max(abs(seg{:,2}))+100];
% [p1,y1]=ksdensity(seg{:,2},'support',[0 range(2)+1e-6]);
% range=[-max(abs(y1)),max(abs(y1))];
% ninterp=200;
% py1=interp1([range(1)-1e-6,y1,range(2)+1e-6],[0,p1,0],linspace(range(1),range(2),ninterp));
% py2=interp1([range(1)-1e-6,-y1,range(2)+1e-6],[0,p1(end:-1:1),0],linspace(range(1),range(2),ninterp));
% py1y2=conv(py1,py2);
% o(1)=sqrt(var(linspace(range(1),range(2),ninterp),px1x2(ninterp:ninterp*2-1))+...
%     var(linspace(range(1),range(2),ninterp),py1y2(ninterp:ninterp*2-1)))/5;
o(1)=mean(mpdm)+3*mean(vpdm);
% mpdm=[];
% for j=4:(3+2*nChannel)
%     for i=2:max(seg{:,1})
%         idx=seg{:,1}==i;
%         idx0=seg{:,1}==i-1;
%         mat=log(seg{idx,j});
%         mat0=log(seg{idx0,j});
%         distmat=pdist2(mat,mat0);
%         num=round(sqrt(length(mat)*length(mat0)));
%         distmat=sort(distmat(:),'ascend');
%         mpdm(i)=max(distmat(1:num));
%         vpdm(i)=std(distmat(1:num));
%         %    resmat=[resmat;[distmat(1:num),zeros(num,1)+i]];
%     end
%     o(j-2)= max(mpdm);   
% end

o(2:(1+2*nChannel))=std(log10(seg{:,4:(3+2*nChannel)}));