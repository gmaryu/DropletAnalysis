function manualMistrackingRemoval(sumTrackTable, dataPath, datadate, position, targets, colNames, mannualEdit, loadSaveFile, isClear)
%% main loop by the position
% name of sub-directory
plotType='cleaned_tracking';
formatOut=30;

for i=1:length(targets)
    posList{1}{i}=['Pos',num2str(targets(i))];
end

posIdx=1;          % for while loop
while posIdx <= length(targets)
    % make new figure object, variable and new folders
    %fig=figure();
    rm_ids=[];
    
    for n=1:size(targets)
        if strcmp(['Pos',num2str(position(posIdx))],posList{1}{n})
            relPos=n; % relative position in genPos
            disp(relPos);
            % make directory for each position to save mat file
            imPath=fullfile(dataPath, plotType, posList{1}{n});
            mkdir(imPath)
        else
        end
    end
    %{
    % load tracking result
    cPos=position(posIdx); % current position
    disp(['You are working on Pos',num2str(cPos)]);
    tmptracks=sumTrackTable{posIdx};
    ori_ids=extractfield(tmptracks,'id'); 
    
    % If 'loadSaveFile' is true, overwrite position index and load the data in progress.
    if posIdx==1 && loadSaveFile==1
        savedFile=load(fullfile(dataPath, plotType,'data_correction_in_progress.mat'),'ii','posIdx','rm_ids');
        %savedFile=load(fullfile(dataPath, plotType,'data_correction_in_progress.mat'),'tmptracks','ii','posIdx','rm_ids');
        %tmptracks=savedFile.tmptracks;
        rm_ids=savedFile.rm_ids;
        posIdx=savedFile.posIdx;
        tmptracks=sumTrackTable{posIdx};
        cPos=position(posIdx); % current position
        disp(['You are working on Pos',num2str(cPos)]);
    end
    
    % loop by droplets
    ii=1;
    while ii <= length(tmptracks)

        if ii==1&&loadSaveFile==1&&~isClear
            % load saved matfile 
            disp(savedFile.ii); % last droplet number in saved mat file
            ii=savedFile.ii; % overwrite 

        else
            if ii==1 && loadSaveFile==1
                ii=savedFile.ii-1;
                
                %disp(ii);
                %disp(['1,load  ',num2str(ii)]);
            end
            
            % set original signal
            eval(['tmpSig=tmptracks(',num2str(ii),').Feat.',colNames{1},';']);

        end
        disp(['ii: ',num2str(ii)]); % droplet id
        %% selection peaks
        plot(tmpSig);
        if mannualEdit==1
            [x,y,button] = ginput(1);
            
            if button=='a'% keep droplet
                disp('ok'); 
            elseif button=='w' % remove droplet
                rm_ids=[rm_ids, ori_ids(ii)]
            elseif button=='s'
                savefile=fullfile(dataPath, plotType, ['data_correction_in_progress_',datestr(now,formatOut),'.mat']);
                save(savefile,'rm_ids','ii','posIdx','-v7.3');
                disp('saving finished');
            elseif button=='q' %exit
                return
            end
        end
        disp([num2str(ii),'/',num2str(size(tmptracks))])
        ii=ii+1;
    end
    ib=ismember(ori_ids,rm_ids);
    tmptracks(ib,:)=[];
    sumTrackTable{posIdx}=tmptracks;
    
    savefile2=fullfile(imPath,'removed_ID.mat');
    save(savefile2,'tmptracks','rm_ids','-v7.3');
    close(fig);
    %}
    posIdx=posIdx+1;
    
end
%%
%save([datadate,'_cleaned_tracks_',datestr(now,formatOut)],'sumTrackTable');