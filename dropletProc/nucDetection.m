function [imnc2,imcy2]=nucDetection(oimg, gB, S1, skwth)

%% memory allocation
imnc=zeros(size(oimg));
imcy=zeros(size(oimg));
imnc2=zeros(size(oimg));
imcy2=zeros(size(oimg));
imd=zeros(size(oimg));
%skwth=[2.0, 1.0, 0.5, 0.3, -0.25];

% draw nuclear detection result
for segid = 1:size(S1) % loop by field number of regoinprep
    tmp1=(oimg==S1(segid).Label); % single droplet in cleaned watershed result
    tmp2=gB.*tmp1; % single droplet in raw image, gray
    pv=S1(segid).PixelValues;
    Th=multithresh(pv,8);

    if S1(segid).skw > skwth(1) % skewness > 2.0
        bwd=imbinarize(tmp2,Th(5));
    elseif skwth(1) > S1(segid).skw && S1(segid).skw > skwth(2) % 2.0 > skewness > 1
        bwd=imbinarize(tmp2,Th(6));
        %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
    elseif skwth(2) > S1(segid).skw && S1(segid).skw > skwth(3) % 1.0 > skewness > 0.5
        bwd=imbinarize(tmp2,Th(7));
        %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
    elseif skwth(3) > S1(segid).skw && S1(segid).skw > skwth(4) % 0.5 > skewness > 0.3
        bwd=imbinarize(tmp2,Th(8));
        se = strel('disk',2);
        bwd=imdilate(bwd, se);
        %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
    elseif skwth(4) > S1(segid).skw && S1(segid).skw > skwth(5) % 0.3 > skewness > -0.25
        bwd=imbinarize(tmp2,Th(8));
        se = strel('disk',2);
        bwd=imopen(bwd, se);
        %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
    else % skewness < -0.25
        bwd=0;
        %disp(segid);
    end
    
    %
    if max(bwd(:)) > 0
        imnc=imnc+bwd; % add nuclear on blank matrix % boolean
    else
        imd=imd+tmp1; % add droplet shape that doesn't contain nuclear on blank matrix
    end
    
    imnc2=oimg.*logical(imnc+imd);
    imcy2=oimg-(oimg.*imnc);
end