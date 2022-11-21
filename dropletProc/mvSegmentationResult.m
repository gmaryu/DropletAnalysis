function mvSegmentationResult(dateList, posList, dataPath, resultPath)
if isempty(dataPath)
    dataPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\data\';
end
if isempty(resultPath)
    resultPath='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\';
end
%% segmentation and intensity quantification
for pp=1:length(dateList)
    for pj=1:length(posList{pp})
        foldername=[dataPath,dateList{pp},'/',posList{pp}{pj}];
        addpath(genpath(pwd))
        cd(foldername)
        
        foldername2=[resultPath,dateList{pp},'/SegmentationResults/',posList{pp}{pj},'_segmentation'];
        if ~isempty(dir('*segment*'))
            movefile('*segment*', foldername2)
        elseif ~isempty(dir('*nuc*'))
            movefile('*_nuc_*',foldername2')
        elseif ~isempty(dir('*cyto*'))
            movefile('*_cyto_*',foldername2')
        end
    end
    
end