function [oscillation_table_nuc, oscillation_table_cyto]=splitOscillationTable(oscillation_table)
%{
    split one oscillation_table to ocillation_table_nuc and
    oscillation_table_cyto
%}
%%
oscillation_table_nuc=oscillation_table;
oscillation_table_cyto=oscillation_table;
for i=1:length(oscillation_table)
    tmp_nuc=oscillation_table_nuc(i).data;
    tmp_cyto=oscillation_table_cyto(i).data;
    
    tmp_nuc(tmp_nuc.Nuc==0,:)=[];
    tmp_cyto(tmp_cyto.Nuc==1,:)=[];
    
    oscillation_table_nuc(i).data=tmp_nuc;
    oscillation_table_cyto(i).data=tmp_cyto;
end

