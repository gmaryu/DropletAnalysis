%function mkDropletMovie(tracks, oscillation_table, input_dropletIDs, position, dataPath, savePath)
%%
%input_dropletIDs=cycB_oscillation_struct; % manually selected droplets and oscillation cycles
oscillation_data=osci_table.data;
input_dropletIDs=unique(oscillation_data.dropID); % automatically detected all oscillations

dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\data\20200127_FRET_NLS_frozen\Pos25';
savePath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20200127_FRET_mRNA\gif\Pos25\CFP\';

% check format
if isstruct(input_dropletIDs)
    droplet_id=input_dropletIDs.dropletID;
elseif isvector(input_dropletIDs)
    droplet_id=input_dropletIDs;
else
    disp('unknown format');
    return;
end

extraframe=10;

%%
nameFormat='img_000000%03d_5-CFP_000.tif';
%nameFormat='img_000000%03d_3-mCherry_000.tif';
%foldername=dataPathTable{relPos};
foldername=dataPath;
%foldername=['\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\data\',dateList{1},'/',posList{1}{relPos}];
%droplet_id=linspace(1,163);

%% add feild for orifinal id info

for segid =1:size(tracks, 1)
    
    if ismember(segid, droplet_id)
        disp(segid);
        xcnt=[tracks(segid).Feat.ycoord];
        ycnt=[tracks(segid).Feat.xcoord];
        tdata=[tracks(segid).Feat.t];
        input_index=find(droplet_id==segid);
        if isstruct(input_dropletIDs)
            
            if iscell(input_dropletIDs.cycleID{input_index})
                cycles=cell2mat(input_dropletIDs.cycleID{input_index});
            else
                cycles=input_dropletIDs.cycleID{input_index};
            end
            
            firstCycle=min(cycles);
            lastCycle=max(cycles);
            tmp=find(oscidrop==droplet_id(input_index)); % table index of target droplet
            tmp_peak=osci.peakTime(tmp); % peak detected timepoints
            
            if firstCycle==1
                firstframe=tracks(droplet_id(input_index)).Feat.t(1);
            else
                firstframe=tmp_peak(firstCycle-1);
            end
            
            if tmp_peak(lastCycle)+extraframe > max(tracks(droplet_id(input_index)).Feat.t)
                lastframe=max(tracks(droplet_id(input_index)).Feat.t);
            else
                lastframe=tmp_peak(lastCycle); %last peak
            end
        elseif isvector(input_dropletIDs)
            firstframe=min(oscillation_data.peakTime(oscillation_data.dropID==segid))+min(tdata);
            lastframe=max(oscillation_data.peakTime(oscillation_data.dropID==segid))+min(tdata);
            %firstframe=min(tracks(droplet_id(input_index)).Feat.t);
            %lastframe=max(tracks(droplet_id(input_index)).Feat.t);

        else
            disp('unknown format');
            return;
        end
        
        h = figure();
        t=firstframe;
        
        if lastframe+extraframe > max(tracks(droplet_id(input_index)).Feat.t)
            endframe=max(tracks(droplet_id(input_index)).Feat.t);
        else
            endframe=lastframe+extraframe; %last peak
        end
        while t < endframe
            positionname=split(dataPath,'\'); positionname=positionname(end); positionname=positionname{1};
            gifname= strcat(positionname,'_TrackingID_',num2str(segid),'.gif');
            X = [positionname, ': TrackID', num2str(segid)];
            fname= sprintf(nameFormat,t);
            Y=[foldername,'/', fname];
            disp(Y);
            %disp(X);
            %X = [X, ' / ', fname];
            
            %disp(gifname);
            
            img= imread(Y);
            gimg= double(img);
            
            % normalize
            ngimg=(gimg-min(gimg(:)))/(max(gimg(:))-min(gimg(:)))*255;
            

            %disp(['t: ', num2str(t)]);
            tindex=find(tdata==t);
            %disp(['time index: ', num2str(tindex)]);
            
            if round(ycnt(tindex))-50 < 1
                tmp=horzcat(zeros(101,abs(round(ycnt(tindex)-50))), ngimg(round(xcnt(tindex)-50):round(xcnt(tindex)+50), 1:round(ycnt(tindex)+50)));
            elseif round(ycnt(tindex))+50 > size(img,1)
                tmp=horzcat(ngimg(round(xcnt(tindex)-50):round(xcnt(tindex)+50), round(ycnt(tindex)-50):end), zeros(101,abs(round(ycnt(tindex)+50-size(img,1)))));
            elseif round(xcnt(tindex))-50 < 1
                tmp=vertcat(zeros(abs(round(xcnt(tindex)-50)),101), ngimg(1:round(xcnt(tindex)+50), round(ycnt(tindex)-50):round(ycnt(tindex)+50)));
            elseif round(xcnt(tindex))+50 > size(img,2)
                tmp=vertcat(ngimg(round(xcnt(tindex)-50):end, round(ycnt(tindex)-50):round(ycnt(tindex)+50)), zeros(abs(round(xcnt(tindex)+50-size(img,1))),101));
            else
                tmp=ngimg(round(xcnt(tindex)-50):round(xcnt(tindex)+50), round(ycnt(tindex)-50):round(ycnt(tindex)+50));
            end

            frame = getframe(h);
            if t == firstframe
                imwrite(tmp,gifname,'gif', 'Loopcount',inf); 
            else
                imwrite(tmp,gifname,'gif','WriteMode','append');
            end
           
            t=t+1;
        end
        close(h);
        %break;
    end
end