function plotColData(plotParam)
%% params
sumTrackTable=plotParam{1};
posList=plotParam{2};
dateList=plotParam{3};
pos=plotParam{4}; % position
sdir=plotParam{5}; % store directory
plotType=plotParam{6}; %plotType
% left Y-axis
lax=plotParam{7}; % colum name for left y-axis
lyl=plotParam{8}; % label for left y-axis
% right Y-axis 
rax=plotParam{9}; % colum name for right y-axis
ryl=plotParam{10};% label for right y-axis
% time
interval=plotParam{11}; % interval
xl=plotParam{12}; % label for x-axis
droplets=plotParam{13}; % droplet IDs ot 'all'
saveExtension=plotParam{14}; % file extension type for save as image
formatOut=30;

%% plotting
for p = 1:size(pos,2)
    % data
    for n = 1:size(posList{1},2)
        if strcmp(['Pos', num2str(pos(p))], posList{1}{n})
            break;
        end
    end
    tracks=sumTrackTable{n};
    spath=fullfile(sdir,dateList{1},plotType{1},[posList{1}{n},'_',datestr(now,formatOut)]);
    disp(spath);
    mkdir(spath);
    cd(spath);

    % plot all droplets data
    if strcmp(droplets,'all')
        for i =1:size(tracks, 1)
            xdata=tracks(i).Feat.t*interval;
            figure('visible','off');hold on;
            title(tracks(i).id);
            %title([tracks(i).DyeSort_ID, tracks(i).DyeMean]);
            %title([tracks(i).id, tracks(i).Feat.BFP_whole(1)]);
            xlabel(xl);
            yyaxis left;
            for c=1:size(lax,2)
                tmpSig1=eval(['tracks(',num2str(i),').Feat.',lax{c}]);
                %disp(lax{c});
                plot(xdata, tmpSig1);
            end
            ylabel(lyl);
            
            if isempty(rax) % only left axis
                name_cells = {strrep(lax, '_',' ')};
            else % both axis
                yyaxis right;
                for d=1:size(rax,2)
                    tmpSig2=eval(['tracks(',num2str(i),').Feat.',rax{d}]);
                    plot(xdata, tmpSig2);
                    name_cells = {strrep(lax, '_',' '),strrep(rax, '_',' ')};
                end
                ylabel(ryl);
            end
            
            % legend
            nameCol = cat(2, name_cells{:});
            %disp(nameCol);
            legend(nameCol);
            
            % save fig as png
            if strcmp(saveExtension,'png')
                fname = ['TrackID_',num2str(tracks(i).id),'.png'];
            elseif strcmp(saveExtension,'eps')
                fname = ['TrackID_',num2str(tracks(i).id),'.eps'];
            else
                disp('This file type is not supported.');
            end
            %fname = ['TrackID_',num2str(tracks(i).DyeSort_ID),'.png'];
            disp(fname);
            saveas(gcf, fname);
            close;
        end
    else
        % plot selected doplets
        for j =1:length(droplets)
            i=droplets(j);
            xdata=tracks(i).Feat.t*interval;
            figure('visible','off');hold on;
            title(tracks(i).id);
            %title([tracks(i).DyeSort_ID, tracks(i).DyeMean]);
            xlabel(xl);
            yyaxis left;
            for c=1:size(lax,2)
                tmpSig1=eval(['tracks(',num2str(i),').Feat.',lax{c}]);
                %disp(lax{c});
                plot(xdata, tmpSig1);
            end
            ylabel(lyl);
            
            if isempty(rax) % only left axis
                name_cells = {strrep(lax, '_',' ')};
            else % both axis
                yyaxis right;
                for d=1:size(rax,2)
                    tmpSig2=eval(['tracks(',num2str(i),').Feat.',rax{d}]);
                    plot(xdata, tmpSig2);
                    name_cells = {strrep(lax, '_',' '),strrep(rax, '_',' ')};
                end
                ylabel(ryl);
            end
            
            % legend
            nameCol = cat(2, name_cells{:});
            %disp(nameCol);
            legend(nameCol);
            
            % save fig as png
            if strcmp(saveExtension,'png')
                fname = ['TrackID_',num2str(tracks(i).id),'.png'];
            elseif strcmp(saveExtension,'eps')
                fname = ['TrackID_',num2str(tracks(i).id),'.eps'];
            else
                disp('This file type is not supported.')
            end
            disp(fname);
            saveas(gcf, fname);
            close;
        end
    end
end