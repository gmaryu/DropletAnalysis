function [ss5, S1]=BFsegmentationDroplets(img,t,resAmp,areaT,EccentricityT,prthresh,cFormat)
% convert to double image
gimg= double(img);

% normalize
ngimg=(gimg-min(gimg(:)))/(max(gimg(:))-min(gimg(:)));

% edge detection for segmentation
c=0.0005;
gimg2=imgaussfilt(ngimg,6*resAmp)-imgaussfilt(ngimg,2*resAmp);
gimg3=gimg2;
gimg3(ngimg<0)=0;

sigma=3*resAmp;
[gDxx,gDxy,gDyy] = Hessian2D(gimg2,sigma);
% Correct for scale
gDxx = (sigma^2)*gDxx;
gDxy = (sigma^2)*gDxy;
gDyy = (sigma^2)*gDyy;

% Calculate (abs sorted) eigenvalues and vectors
[gLambda2,gLambda1,~,~]=eig2image(gDxx,gDxy,gDyy);
gS2 = gLambda1.^2;% + Lambda2.^2;
gI2=(ones(size(gimg2))-exp(-gS2/c)).*(gLambda1<0).*gimg3;
tmpthresh=graythresh(gI2)/4;
ngimg=gI2>tmpthresh;
ngimgl=bwlabel(ngimg); %binalize and labling
%figure; imshow(ngimgl);
%saveas(gcf,'ngimgl.tif');

gS=regionprops(ngimgl,'Area','Eccentricity');
delidx=([gS.Area]<areaT(1)&[gS.Eccentricity]<EccentricityT)|([gS.Area]<areaT(2)&[gS.Eccentricity]<EccentricityT); % pick-up of junk labels

% replace labels to zero
for iij=1:length(delidx)
    if delidx(iij)
        ngimgl(ngimgl==iij)=0;
    end
end
%figure; imshow(ngimgl);
ngimg=logical(ngimgl); % binalize of label image
ngimg=gimg3&ngimg;
seeds=imextendedmax(bwdist(ngimg),2);
gimg4=max(gimg(:))- gimg;
gss1=watershed(imimposemin(gimg4,seeds));
gss2=watershed(imimposemin(gimg4,gss1==0|seeds));
%figure; imagesc(gss2);
gss2(gss2==1)=0;
gss3=logical(gss2);
gss3=imopen(gss3,strel('disk',round(2*resAmp)));
gss2=double(gss2).*double(gss3);
SE = strel('disk',2);
gss2=imerode(gss2,SE);

fname= sprintf(cFormat{1},t);
B=imread(fname);
gB=double(B);
S1=regionprops(gss2,B,'Centroid','Area','Perimeter',...
    'PixelIdxList','PixelValues','Eccentricity','EulerNumber','MeanIntensity');

% collection and overwrite of junk label
%disp(size(S1));
delid=[S1.Area]>areaT(2)|[S1.Area]<areaT(1); % too big or too small
delid=delid|[S1.Eccentricity]>EccentricityT; % roundness
delid=delid|[S1.Perimeter].^2./[S1.Area]>(1+prthresh)*4*pi; % mismutch of perimeter and area ?
delid=delid|[S1.EulerNumber]>1|[S1.EulerNumber]<1;
%disp(sum(delid));
S1(delid==1)=[];

% apply delid to watershed result (cleaned segmentation result)
%figure; imagesc(gss2);
ss5=gss2;
for ij=1:max(ss5(:))
    %for ij=2:max(gss2(:))
    if delid(ij)
        ss5(ss5==ij)=0;
    end
end

ss5l=unique(ss5);
uqss5=num2cell(ss5l(2:end));
[S1.Label]=uqss5{:};
end