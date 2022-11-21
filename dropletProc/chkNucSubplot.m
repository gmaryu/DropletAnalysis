function chkNucSubplot(i,sg,mkimg,S1,smplist,cvth1,cvth2,stdth,skwth,savedir)
%%
for did = 1:size(S1) % loop by field number of regoinprep
    for n=1:length(smplist)
        if did == smplist(n)
            %%
            segid=did;
            pv=S1(segid).PixelValues;
            meanInt=S1(segid).MeanIntensity;
            skw=S1(segid).skewness; % default:0
            stdInt=S1(segid).stdInt;
            cv=stdInt/meanInt;            
            Th=multithresh(pv,7); % multiple thresh
            tmpl=(sg==S1(segid).Label); % single droplet in cleaned watershed result (binary)
            tmp2=mkimg.*tmpl; % single droplt iamge with raw pixel value
            cropsize=40;
            
            if cv > cvth1 && stdInt > stdth
                if skw > skwth
                    hfig=figure('Position', [10 10 1800 600]);
                    sgtitle(['StdInt: ',num2str(stdInt), '/ segid:', num2str(segid),'/ CV:',num2str(cv)]);
                    % imscale droplet raw image
                    subplot(2,5,1);
                    cnt=S1(segid).Centroid;
                    if round(cnt(1)-cropsize) < 0
                        tmp=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        tmp=horzcat(tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        tmp=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), tmp2(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        tmp=vertcat(mkimg(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        tmp=tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    imagesc(tmp);
                    set(gca,'dataAspectRatio',[1 1 1])
                    
                    % draw intensity distribution histogram
                    subplot(2,5,6);
                    histogram(double(pv));
                    title(skw);
                    
                    disp('Clear nuc');
                    % draw overlay image
                    subplot(2,5,2);
                    bwd=imbinarize(tmp2,Th(1));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,3);
                    bwd=imbinarize(tmp2,Th(2));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,4);
                    bwd=imbinarize(tmp2,Th(3));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,5);
                    bwd=imbinarize(tmp2,Th(4));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,7);
                    bwd=imbinarize(tmp2,Th(5));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,8);
                    bwd=imbinarize(tmp2,Th(6));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,9);
                    bwd=imbinarize(tmp2,Th(7));
                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);
                    
                    subplot(2,5,10);
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

                    if round(cnt(1)-cropsize) < 0
                        bwd2=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        bwd2=horzcat(bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        bwd2=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), bwd(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        bwd2=vertcat(bwd(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        bwd2=bwd(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    olimg=labeloverlay(tmp3, bwd2);
                    imshow(olimg);

                else
                    hfig=figure();
                    sgtitle(['StdInt: ',num2str(stdInt), '/ segid:', num2str(segid),'/ CV:',num2str(cv)]);
                    % imscale droplet raw image
                    subplot(1,2,1);
                    cnt=S1(segid).Centroid;
                    if round(cnt(1)-cropsize) < 0
                        tmp=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                    elseif round(cnt(1)+cropsize) > size(sg,1)
                        tmp=horzcat(tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                    elseif round(cnt(2)-cropsize) < 0
                        tmp=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), tmp2(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                    elseif round(cnt(2)+cropsize) > size(sg,2)
                        tmp=vertcat(tmp2(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                    else
                        tmp=tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                    end
                    tmp3=uint8((tmp / max(tmp(:)) * 255));
                    imshow(tmp3);
                    
                    % draw intensity distribution histogram
                    subplot(1,2,2);
                    histogram(double(pv));
                    title(skw);
                    
                end
                    

            else
                %disp('Dim nuc');
                hfig=figure();
                sgtitle(['StdInt: ',num2str(stdInt), '/ segid:', num2str(segid),'/ CV:',num2str(cv)]);
                % imscale droplet raw image
                subplot(1,2,1);
                cnt=S1(segid).Centroid;
                if round(cnt(1)-cropsize) < 0
                    tmp=horzcat(zeros(cropsize*2+1,abs(round(cnt(1)-cropsize))), tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), 1:round(cnt(1)+cropsize)));
                elseif round(cnt(1)+cropsize) > size(sg,1)
                    tmp=horzcat(tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):end), zeros(cropsize*2+1,abs(round(cnt(1)+cropsize-size(sg,1)))));
                elseif round(cnt(2)-cropsize) < 0
                    tmp=vertcat(zeros(abs(round(cnt(2)-cropsize)),cropsize*2+1), tmp2(1:round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize)));
                elseif round(cnt(2)+cropsize) > size(sg,2)
                    tmp=vertcat(tmp2(round(cnt(2)-cropsize):end, round(cnt(1)-cropsize):round(cnt(1)+cropsize)), zeros(abs(round(cnt(2)+cropsize-size(sg,1))),cropsize*2+1));
                else
                    tmp=tmp2(round(cnt(2)-cropsize):round(cnt(2)+cropsize), round(cnt(1)-cropsize):round(cnt(1)+cropsize));
                end
                tmp3=uint8((tmp / max(tmp(:)) * 255));
                imshow(tmp3);
                
                % draw intensity distribution histogram
                subplot(1,2,2);
                histogram(double(pv));
                title(skw);
                
            end
            nameFormat1='subplot3_000000%03d_';
            fname2= sprintf(nameFormat1,i);
            %fn=[savedir,fname2,num2str(S1(did).Label),'.png'];
            fn=[savedir,fname2,num2str(segid),'.png'];
            if ~exist(savedir, 'dir')
                mkdir(savedir)
            end
            %disp(fn);
            %fn=['\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20200214_pGM020-pGM021-pGM023\Pos18_nuc2\',fname2,num2str(S1(segid).Label),'.png'];
            %print(hfig,'-dpng','-r0',fn);
            saveas(gcf,fn);
            close(hfig)
        end
    end
    
end
end