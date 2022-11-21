function [tracks,mean_ctv]=linearBaselineCorrection8(tracks, osci_table,colNames,tentative_trough_frames)
%% linear correction of baseline
%{
Description of function
    Estimate linear baseline with trough values for no-nucleus droplets.
    Corrected signal values are calculated by adding estimated baseline values to original signal values in each droplet.
    Corrected signal values are also normalized with a mean value of 1st trough value of all no-nucleus droplets.
tracks: tracking result object for target position (e.g. sumtrackTable{1})
osci_table: peak detection result
colNames: list of column names 
colNames{1}:target column e.g. *_ratio_whole
%}
    
    %% data correction by estimated value
    if length(colNames)==1 
        nonuc_1st_trough=[];
        for n =1:size(tracks)
            %%
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
      
            if length(trough_v) >= 2 % fitting with original data
                trough_tp1=trough_tp(1);
                peak_tp1=peak_tp(1);
                peak_tp2=peak_tp(2);
                % linear fitting
                lm=fitlm(tbl, 'linear');
                est_b=lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*tp; % estimated baseline
                est_icpt=lm.Coefficients.Estimate(1)*ones(length(tp),1); % estimated intercept vector
                est_res=est_icpt-est_b; % residue (difference)
                est_val=tmpSig+est_res; % estimated value

            elseif length(trough_v) == 1 %
                trough_tp1=trough_tp(1);
                peak_tp1=peak_tp(1);
                peak_tp2=peak_tp(2);
                trough_tp2=peak_tp2+tentative_trough_frames+min(tp)-1; %  
                if trough_tp2 > length(tmpSig)
                    trough_tp2 = length(tmpSig);
                end
                trough_tp2_v=tmpSig(trough_tp2);
                tbl.trough_tp(end+1)=trough_tp2;
                tbl.trough_v(end)=trough_tp2_v;
                
                % linear fitting
                lm=fitlm(tbl, 'linear');
                est_b=lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*tp; % estimated baseline
                est_icpt=lm.Coefficients.Estimate(1)*ones(length(tp),1); % estimated intercept vector
                est_res=est_icpt-est_b; % residue (difference)
                est_val=tmpSig+est_res; % estimated value               
            else
                est_val=NaN*ones(length(tp),1);
                est_val2=NaN*ones(length(tp),1);
                est_b=NaN*ones(length(tp),1);
            end
            %%
            
            min_v_1st_cycle=min(est_val(peak_tp1:peak_tp2));
            nonuc_1st_trough=[nonuc_1st_trough,min_v_1st_cycle];
            %%
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear71=est_val;']);
            tracks(n).Feat.est_linear_baseline=est_b;
        end
        
        %% set average trough value to 1.0
        %disp(corrected_trough_values);
        mean_ctv=nanmean(nonuc_1st_trough);

        for n =1:size(tracks)
            % get droplet and estimated data
            eval(['did=tracks(',num2str(n),').id;']);
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},'_linear71;']);
            eval(['tp=tracks(',num2str(n),').Feat.t;']);
            %disp(n);
            %disp(size(mean_estb));
            tmpSig=tmpSig./mean_ctv;
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear72=tmpSig;']);
        end
        
    else
        disp("You need to select one target column as colNames.");
    end
    
end

%}