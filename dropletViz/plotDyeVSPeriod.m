function plotDyeVSPeriod(position, oscillation_table)
%{
    plot a scatter for "DyeIntensity" vs "Period length" in each cycle (not include 1st cycle)
    position: a list of positions where to analyze
    oscillation_table:
%}
%% variables to make a plot
PosX_period=[];
PosX_cycleID=[];
PosX_dye=[];

%%
allpositions=extractfield(oscillation_table,'position');
selected_position=ismember(allpositions, position);

%% main
for p =1:length(oscillation_table)
    if selected_position(p)==1
        % select temporal position's information
        PosN_osci=oscillation_table(p).data;
        
        % select features
        periodpeak=PosN_osci.periodPeak;
        peaktime=PosN_osci.peakTime;
        cycleID=PosN_osci.cycleID;
        dropletID=PosN_osci.dropID;
        dyeInt=PosN_osci.DyeIntensity;
        
        for j=1:length(dropletID)
            % cycle info
            did=dropletID(j);
            cid=cycleID(j);
            pp=periodpeak(j);
            pt=peaktime(j);
            di=dyeInt(j);
            
            % make vector for plot
            PosX_cycleID=vertcat(PosX_cycleID, cid);
            PosX_dye=vertcat(PosX_dye, di);
            if isnan(pp) && cid==1
                PosX_period=vertcat(PosX_period, pt);
            else
                PosX_period=vertcat(PosX_period, pp);
            end
        end
        
    end
end
figure();
scatter(PosX_dye, PosX_period, 5, PosX_cycleID,'filled');
colorbar;
xlabel('Dye Intensity');
ylabel('Period length (Frame)');