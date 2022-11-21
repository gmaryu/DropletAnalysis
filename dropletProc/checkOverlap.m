function [out,ov]=checkOverlap(in,c,r,n,overthresh)
thresh=overthresh;
area=r*r*pi;
count=0;
out=in;
for i=max(1,round(c(1)-r)):min(size(in,1),round(c(1)+r))
    for j=max(1,round(c(2)-r)):min(size(in,2),round(c(2)+r))
        if (i-c(1))*(i-c(1))+(j-c(2))*(j-c(2))<r*r
           if out(i,j)==0
               out(i,j)=n;
           else
               count=count+1;
               if count>area*thresh
                   out=in;
                   ov=1;
                   return;
               end
           end
        end
    end
end
ov=0;
        