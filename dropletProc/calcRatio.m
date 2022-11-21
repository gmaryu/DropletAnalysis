function tracks=calcRatio(tracks,numeratorCh,denominatorCh,colName)

for i=1:size(tracks)
    eval(['ch1=tracks(',num2str(i),').Feat.',numeratorCh,';']);
    eval(['ch2=tracks(',num2str(i),').Feat.',denominatorCh,';']);
    eval(['tracks(',num2str(i),').Feat.',colName,'=ch1./ch2;']);
end
end