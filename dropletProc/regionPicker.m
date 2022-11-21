function sroi=regionPicker(datafolder, dateList, posList, nameFormat, frameNum)
%%%%
% selection of a ROI for background subtraction
% minimum projection for ROI selection
%%%%
%%
nf=100;
for p=1:size(posList,2)
    filename1=sprintf(nameFormat,0);
    fileLoc1=[datafolder, dateList,'/',posList{p},'\',filename1];
    temp=imread(fileLoc1);
    imgsize=size(temp);
    imgstk=zeros(imgsize(1), imgsize(2), nf);
    tpvec=linspace(1,frameNum,nf);
    for q=1:numel(tpvec)
        disp(q);
        fname= sprintf(nameFormat,round(tpvec(q)));
        fileLoc=[datafolder, dateList,'/',posList{p},'\',fname];
        tempimg=imread(fileLoc);
        imgstk(:,:,q)=tempimg;
    end
    minplane=min(imgstk,[],3);
    %disp(minplane);
    f(p)=figure(p);
    %temp=mean(cat(3,double(temp1),double(temp2),double(temp3)),3);
    imshow(minplane,[min(minplane(:)), prctile(minplane(:),95)]);
    title(posList{p});
    roi=drawrectangle;
    sroi{p}=roi.Vertices;
    close(f(p));
end