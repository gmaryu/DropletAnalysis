function [img2, meanbg]=bgSubtraction(img, roicorr)
meanbg=mean(img(int64(roicorr(1,1)):int64(roicorr(3,1)), int64(roicorr(1,2)):int64(roicorr(2,2))),'all');
img2=img-meanbg;
img2(img2<0)=0;
%disp(min(min(img2)));
end
