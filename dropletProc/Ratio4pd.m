function tracks=Ratio4pd(tracks,colName,mag)

for i=1:size(tracks)
    eval(['ch1=tracks(',num2str(i),').Feat.',colName,';']);
    eval(['tracks(',num2str(i),').Feat.sig4pd=ch1*mag;']);
end
end