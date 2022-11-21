function [oscillation_table_nuc, oscillation_table_cyto]=classifyOsciByNuc(oscillation_table)
%{
Description
    split oscillation_table into two groups (nuc or cyto)
 
Args
    oscillation_table: genarated table by peak detection function
%}
%%
% new tabels for return value
oscillation_table_nuc=oscillation_table;
oscillation_table_cyto=oscillation_table;

for i =1:length(oscillation_table)
    oscis_nuc=oscillation_table(i).data;
    oscis_nuc(oscis_nuc.Nuc == 0, :) = [];
    oscis_cyto=oscillation_table(i).data;
    oscis_cyto(oscis_cyto.Nuc == 1, :) = [];
    oscillation_table_nuc(i).data=oscis_nuc;
    oscillation_table_cyto(i).data=oscis_cyto;
end

end 