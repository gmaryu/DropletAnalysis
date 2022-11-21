function tracks=classifyNucCyt(tracks)
%{
    Description
        Detect the area difference between nucleus area and cytoplasm area.
        Retuern the frame number that nucleus reagion was detected. 
%}

for did=1:size(tracks)
   disp(did);
   area_nuc=tracks(did).Feat.Nuc_area;
   area_cyt=tracks(did).Feat.Cyto_area;
   cmp_res=area_nuc < area_cyt;
   tracks(did).Feat.Nuc=cmp_res;
   tracks(did).Nuc=sum(cmp_res);

end

end