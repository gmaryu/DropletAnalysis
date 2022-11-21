function tracks=linearBaselineCorrection2(tracks, osci_table, colNames)
%% linear correction of baseline
%{
Description of function
    Estimate linear baseline with trough values. Corrected trough values
    are set around 1.0.
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
            
            if length(trough_v) >= 2
                % linear fitting
                lm=fitlm(tbl, 'linear');
                est_b=lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*tp; % estimated baseline
                est_icpt=lm.Coefficients.Estimate(1)*ones(length(tp),1); % estimated intercept vector
                est_res=ones(length(est_b),1)-est_b; % residue (difference)
                est_val=tmpSig+est_res; % estimated value
               
            else
                est_val=NaN*ones(length(tp),1);
                est_b=NaN*ones(length(tp),1);
            end
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear2=est_val;']);
            tracks(n).Feat.est_linear_baseline=est_b;
        end
    end
end

%}
