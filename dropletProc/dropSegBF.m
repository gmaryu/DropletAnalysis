function [oimg, S1]=dropSegBF(i, fname, cFormat, segParam)
sigma=segParam{1}; % sigma for droplet segmentation
resAmp=segParam{2}; % resAmp, resolution
minDopDIa=segParam{3}; % minDopDia= % minimum drop diameter?
EccentricityT=segParam{4}; %EccentricityT
prthresh=segParam{5}; %prthresh=0.35; % perimeter thresh?
lenThresh=segParam{6}; %lenThresh=270; % minimum tracking length
areaT=segParam{7}; % areaT= % min and max threshold of area size
areaChangePercent=segParam{8}; % how much area change is allowed
formatOut=30; % date format for 'yyyymmddTHHMMSS'

% load bright field image
img= imread(fname);
gimg= double(img);

% normalize
ngimg=(gimg-min(gimg(:)))/(max(gimg(:))-min(gimg(:)));

% edge detection for segmentation
c=0.0005;
gimg2=imgaussfilt(ngimg,6)-imgaussfilt(ngimg,2);
gimg3=gimg2;
gimg3(ngimg<0)=0;

%sigma=3;
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
delidx=[gS.Area]<(50*resAmp^2)|(([gS.Area]<1000*resAmp^2)&[gS.Eccentricity]<0.95); % pick-up of junk labels

% replace labels to zero
for iij=1:length(delidx)
    if delidx(iij)
        ngimgl(ngimgl==iij)=0;
    end
end
%figure; imshow(ngimgl);
ngimg=logical(ngimgl); % binalize of label image
gimg2=imclose(ngimg,strel('disk',round(5*resAmp))); % morphology filter for binalized image
gimgo=bwareaopen(gimg2,1500); % delete small objects
%{bainari dakede segmentation dekirunodeha?}
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

fname= sprintf(cFormat{1},i-1);
B=imread(fname);
gB=double(B);
S1=regionprops(gss2,B,'Centroid','Area','Perimeter',...
    'PixelIdxList','PixelValues','Eccentricity','EulerNumber');

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
oimg=gss2;
for ij=1:max(oimg(:))
    %for ij=2:max(gss2(:))
    if delid(ij)
        oimg(oimg==ij)=0;
    end
end

ss5l=unique(oimg);
uqss5=num2cell(ss5l(2:end));
[S1.Label]=uqss5{:};
%%
%{
% test code for label and regionprops is coupled.
smplist=randi([1,size(S1,1)],1,10);
disp(smplist);
for n=1:length(smplist)
    tar=S1(smplist(n)).Label;
    cnt=S1(smplist(n)).Centroid;
    disp(cnt);
    figure; imagesc(ss5==tar);
    hold on;
    scatter(cnt(1), cnt(2), 150,'d');
end
%}
stdInt=zeros(size(S1,1),1); % standard dev of pixel value
skw=zeros(size(S1,1),1); % skewness of pixel value
for jj=1:size(S1) % why start from 2
    %for jj=2:size(S1) % why start from 2
    S1(jj).stdInt=std(double(S1(jj).PixelValues));
    S1(jj).skw=skewness(double(S1(jj).PixelValues));
end

end