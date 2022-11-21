function osci_table=readOscillationData(dataPath, pos)
%{
 read oscillation data (mat) that is peakSelection result.
 rename osci -> PosX_osci for organize several position data
%}
osci_table=struct;
posIdx=1;
while posIdx <= length(pos)
    matPath=[dataPath,'\Pos', int2str(pos(posIdx))];
    disp(matPath);
    saveFileName=ls([matPath,'\*.mat']);
    load(fullfile(matPath,saveFileName));
    osci_table(posIdx).position=pos(posIdx);
    osci_table(posIdx).data=osci;
    %eval(['Pos',int2str(pos(posIdx)),'_osci=osci']);
    %clear('osci')
    posIdx=posIdx+1;
end
end