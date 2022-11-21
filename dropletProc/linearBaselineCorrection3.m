function tracks=linearBaselineCorrection3(tracks, osci_table, colNames)
%% linear correction of baseline
%{
Description of function
    Estimate linear baseline with trough values. 
    Corrected signal values are calculated by adding estimated baseline values to original signal values.
    Corrected signal values are also normalized with trough value of 1st value.
tracks: tracking result object for target position (e.g. sumtrackTable{1})
osci_table: peak detection result
colNames: list of column names 
colNames{1}:target column e.g. *_ratio_whole
%}
    
    
    %% data correction by estimated value
    if length(colNames)==1    
        for n =1:size(tracks)
            % get droplet data
            eval(['did=tracks(',num2str(n),').id;']);
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},';']);
            eval(['tp=tracks(',num2str(n),').Feat.t;']);
            
            % get trough data
            troug_idx=find(osci_table.data.dropID==did);
            trough_v=osci_table.data.troughValue(troug_idx);
            trough_v=trough_v(~isnan(trough_v));
            trough_tp=osci_table.data.troughTime(troug_idx);
            trough_tp=trough_tp(~isnan(trough_tp));
            tbl=table(trough_tp,trough_v);
            
            % get peak data
            peak_idx=find(osci_table.data.dropID==did);
            peak_v=osci_table.data.peakValue(peak_idx);
            peak_v=peak_v(~isnan(peak_v));
            peak_tp=osci_table.data.peakTime(peak_idx);
            peak_tp=peak_tp(~isnan(peak_tp));

            
            if length(trough_v) >= 2
                trough_tp1=trough_tp(1);
                peak_tp1=peak_tp(1);
                peak_tp2=peak_tp(2);
                % linear fitting
                lm=fitlm(tbl, 'linear');
                est_b=lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*tp; % estimated baseline
                est_icpt=lm.Coefficients.Estimate(1)*ones(length(tp),1); % estimated intercept vector
                est_res=est_icpt-est_b; % residue (difference)
                est_val=tmpSig+est_res; % estimated value
                min_v_1st_cycle=min(est_val(peak_tp1:peak_tp2));
                est_val2=est_val./min_v_1st_cycle;
               
            else
                est_val=NaN*ones(length(tp),1);
                est_val2=NaN*ones(length(tp),1);
                est_b=NaN*ones(length(tp),1);
            end
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear=est_val;']);
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear3=est_val2;']);
            tracks(n).Feat.est_linear_baseline=est_b;
        end
    end
end

%}