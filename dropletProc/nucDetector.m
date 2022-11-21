function [imnc,imcy,imnc2,imcy2]=nucDetector(sg,mkimg,S1,cvth1,cvth2,stdth,skwth)
imnc=zeros(size(sg));
imcy=zeros(size(sg));
imnc2=zeros(size(sg));
imcy2=zeros(size(sg));
imd=zeros(size(sg));

%%
for segid = 1:size(S1) % loop by field number of regoinprep
    pv=S1(segid).PixelValues;
    meanInt=S1(segid).MeanIntensity;
    skw=S1(segid).skewness;
    stdInt=S1(segid).stdInt;
    cv=stdInt/meanInt;
    Th=multithresh(pv,7); % multiple thresh
    tmpl=(sg==S1(segid).Label); % single droplet in cleaned watershed result (binary)
    reg_tmpl=regionprops(tmpl);
    %disp(reg_tmpl)
    tmp2=mkimg.*tmpl; % single droplt iamge with raw pixel value
    bwd=zeros(size(sg));
    if cv > cvth1 && stdInt > stdth && length(reg_tmpl)==1
        if skw > skwth
            %% nuclear region selection
            if cv > cvth2(1) % cv > 1.5
                bwd=imbinarize(tmp2,Th(1));
            elseif cvth2(1) > cv && cv > cvth2(2) % 2.0 > skewness > 1
                bwd=imbinarize(tmp2,Th(2));
                %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
            elseif cvth2(2) > cv && cv > cvth2(3) % 1.0 > skewness > 0.5
                bwd=imbinarize(tmp2,Th(3));
                %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
            elseif cvth2(3) > cv && cv > cvth2(4) % 0.5 > skewness > 0.3
                bwd=imbinarize(tmp2,Th(4));
                %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
            elseif cvth2(4) > cv && cv > cvth2(5) % 0.3 > skewness > -0.25
                bwd=imbinarize(tmp2,Th(5));
                %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
            elseif cvth2(5) > cv && cv > cvth2(6) % 0.3 > skewness > -0.25
                bwd=imbinarize(tmp2,Th(6));
                %imnc=imnc+bwd; % overlay nuclear position on blank matrix % boolean
            else % skewness < -0.25
                bwd=imbinarize(tmp2,Th(7));
                %disp(segid);
            end
        end
    else
        %disp('dim nuc');
    end 

    %disp(max(bwd(:)));
    if max(bwd(:)) > 0
        reg_bwd=regionprops(bwd);
        reg_bwd_areamax=max(extractfield(reg_bwd,'Area'));
        %disp(reg_bwd)
        if reg_bwd_areamax < reg_tmpl.Area*0.98
            imnc=imnc+bwd; % add nuclear on blank matrix % boolean
        else
            imd=imd+tmpl;
        end
    else
        %disp(segid);
        imd=imd+tmpl; % add droplet shape that doesn't contain nuclear on blank matrix
    end
    
    imnc2=sg.*logical(imnc+imd);
    imcy2=sg-(sg.*imnc);
end
end