function p=calcFeatProb(feat1,feat2,type,par,idx)
% lzd 06-17-2015
% main probability calculation function
% when moving feat1->t-1; feat2->t
% when divide of seperation feat1->t-1 feat2,feat3->t
% when overlap feat1,feat2->t-1, feat3->t
% type is the event string
% par is the std paramter,
% idx is track state(currently is 1)
% numC is number of channel
% return single number p
if idx
    pMove=par(4)*2;
    pLost=eps;
    pMoveIn=par(8);
    pMoveOut=par(9);
    pDisAppear=par(10);
    pAppear=par(11);
    
else
end
switch type
    case 'move'
        if idx
            dist=sqrt((feat2(:,2)-feat1(:,2))^2+(feat2(:,3)-feat1(:,3))^2);
            %diffMat=[dist,abs(log10(feat2(4))-log10(feat1(4))),abs(feat2(5)-feat1(5)),abs(feat2(6)-feat1(6))];
            %p=computeP(diffMat,par(1:4))*pMove+eps;
            diffMat=[dist,abs(log10(feat2(4))-log10(feat1(4))),abs(feat2(5)-feat1(5))];
            p=computeP(diffMat,par(1:3))*pMove+eps;
        else
        end
    case 'moveout'
        if idx
            p=computeP(feat1*4,par(1))+eps+pMoveOut;
        else
        end
        
    case 'movein'
        if idx
            p=computeP(feat1*4,par(1))+eps+pMoveIn;            
            %p=(exp(-5*feat1.^2/par(1).^2)*pBound)+eps+pMoveIn;
        else
        end
        
    case 'disappear'
        if idx
            p=pLost+pDisAppear;
        else
        end
        
    case 'appear'
        if idx
            p=pLost+pAppear;
        else
        end
        
end
%p=log10(p);