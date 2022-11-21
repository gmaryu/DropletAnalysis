function dataTransfer(mode, data, sdir, genPos, ddir)
%{
Description
    move data 
Args
    mode: 'cp' or 'mv'
    data: cell array for part of filename (e.g. {'*nuc*','cyto'}, {'*GFP_000.tif'})
    sdir: path to source data folder 
        (e.g. '\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\data\20190830_CycB1-YFP_Minjun\') 
    genPos: list of the number of position
        (e.g. [0,1,2,3,4,5,6,7])
    ddir: path to destination data folder
        (e.g. '\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\20190830_CycB1-YFP_Minjun\segmantation_results\')
%}
for i=1:length(genPos)
    posList{1}{i}=['Pos',num2str(genPos(i))];
    dpath=fullfile(ddir,posList{1}{i});
    mkdir(dpath);
    
    foldername=[sdir,posList{1}{i}];
    %addpath(genpath(pwd))
    %cd(foldername)
    for n=1:length(data)
        if strcmp(mode,'cp')
            %disp('copying data')
            copyfile(fullfile(foldername,data{n}), dpath);
        elseif strcmp(mode,'mv')
            %disp('move data')
            disp(fullfile(foldername,data{n}));
            disp(dpath);
            movefile(fullfile(foldername,data{n}), dpath);
        end
    end
end
end
