
for pos = 1:length(dataPathTable)
    tmpPath=dataPathTable{pos};
    tmpPos=posList{1}(pos);
    
    gifname= strcat(datadate,'_',tmpPos,'_tracking_',datestr(now,formatOut),'.gif');
    %gifname='aaa';
    disp(gifname);
    
    nf= 400; %get the number of timepoint
    %nameFormat='';
    cd(tmpPath);
    %% Plot tracking result from .mat file
    
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    for i = 1:nf
        fname= sprintf(nameFormat,i-1);
        A=imread(fname);
        %A=double(A);
        %A=(A-min(A(:)))/(max(A(:))-min(A(:)));
        imshow(A);
        
        xpos_g=sumResultTable{pos}.xcoord(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray>0); % good
        ypos_g=sumResultTable{pos}.ycoord(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray>0); % good
        id_g=sumResultTable{pos}.track(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray>0);
        text(xpos_g, ypos_g, cellstr(string(id_g)),'FontSize',6,'FontWeight','bold','Color','g') % good tracking result
        
        xpos_b=sumResultTable{pos}.xcoord(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray==0); % bad
        ypos_b=sumResultTable{pos}.ycoord(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray==0); % bad
        id_b=sumResultTable{pos}.track(sumResultTable{pos}.t==i&sumResultTable{pos}.delArray==0);
        text(xpos_b, ypos_b, cellstr(string(id_b)),'FontSize',6,'FontWeight','bold','Color','r');
        drawnow
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if i == 1
            imwrite(imind,cm,gifname{1},'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,gifname{1},'gif','WriteMode','append');
        end
    end
end
%%

% figure for trajectory
%{
            figure(2)
            hold on
            if i~=1
                tidx=find(sumResultTable{1}.t==i);
                for jj=1:length(tidx)
                    trackIdx=trackArray(tidx(jj));
                    p=find(trackArray==trackIdx&sumResultTable{1}.t==i-1);
                    if ~isempty(p)
                        plot([sumResultTable{1}.xcoord(p),sumResultTable{1}.xcoord(tidx(jj))],...
                            [sumResultTable{1}.ycoord(p),sumResultTable{1}.ycoord(tidx(jj))],'Color',cmap(i,:))
                        hold on
                    end
                end
            end
            drawnow
            frame1 = getframe(gg1);
%}

% figure blended above 2 figures
%{
            gg2=figure(3);
            imshow(A,[])
            frame2 = getframe(gg2);
            imshowpair(frame2im(frame1),frame2im(frame2),'blend');
            frame = getframe(gg2);
            if i>1
                %im2(:,:,:,i) = frame2im(frame);
                %im3(:,:,:,i) = frame2im(frame1);
            else
                temp=frame2im(frame);
                %im2=zeros([size(temp),nf],'uint8');
                %im2(:,:,:,1)=frame2im(frame);
                temp=frame2im(frame1);
                %im3=zeros([size(temp),nf],'uint8');
                %im3(:,:,:,1)=frame2im(frame1);
                
            end
            
        end
%}
%{
        sumResultTable{1}.track=delArray;
        cd ..
        close(h)
        %close(gg1)
        %close(gg2)
        v = VideoWriter([posList{pp}{pj},'_tracking_',datestr(now,formatOut),'.avi']);
        v.FrameRate=6;
        open(v)
        writeVideo(v,im)
        close(v)
%}
%         v = VideoWriter([posList{pp}{pj},'dragonTail_',datestr(now,formatOut),'.avi']);
%         v.FrameRate=6;
%         open(v)
%         writeVideo(v,im2)
%         close(v)
%
%         v = VideoWriter([posList{pp}{pj},'overlay_',datestr(now,formatOut),'.avi']);
%         v.FrameRate=6;
%         open(v)
%         writeVideo(v,im3)
%         close(v)