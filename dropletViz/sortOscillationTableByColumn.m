function oscillation_table=sortOscillationTableByColumn(oscillation_table,sortColumnName)
for t=1:size(oscillation_table,2)
    tempTable=oscillation_table(t).data;
    tempTable=sortrows(tempTable,sortColumnName);
    oscillation_table(t).data=tempTable;
end