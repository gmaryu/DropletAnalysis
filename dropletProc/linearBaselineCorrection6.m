function [tracks,corrected_baseline]=linearBaselineCorrection6(tracks, osci_table,frameNum,colNames)
%% linear correction of baseline
%{
Description of function
    Estimate linear baseline with trough values. 
    Corrected signal values are calculated by adding estimated baseline values to original signal values.
    Corrected signal values are also normalized with a mean value of trough value of all cyctoplasmic droplets.
tracks: tracking result object for target position (e.g. sumtrackTable{1})
osci_table: peak detection result
colNames: list of column names 
colNames{1}:target column e.g. *_ratio_whole
%}
    
    
    %% data correction by estimated value
    if length(colNames)==1    
        corrected_trough_values=[];
        corrected_baseline=[];
        max_data_len=max(extractfield(tracks,'Len'));
        for n =1:size(tracks)
            % get droplet data
            
            eval(['did=tracks(',num2str(n),').id;']);
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},';']);
            eval(['tp=tracks(',num2str(n),').Feat.t;']);
            eval(['nuc=tracks(',num2str(n),').Feat.Nuc;']);
            
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
            
            if nuc < 10
                if length(trough_v) >= 2
                    % linear fitting
                    lm=fitlm(tbl, 'linear');
                    est_b=lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*tp; % estimated baseline
                    est_icpt=lm.Coefficients.Estimate(1)*ones(length(tp),1); % estimated intercept vector
                    est_res=est_icpt-est_b; % residue (difference)
                    est_val=tmpSig+est_res; % estimated value
                    est_b2=ones(frameNum,1)*NaN;
                    est_b2(tp(1):tp(end))=est_b;
                    corrected_baseline=horzcat(corrected_baseline,est_b2);
                    % if current droplet is no-nuc droplet, get trough value
                    % for data correction
                    
                    for p=1:length(peak_tp)-1
                        p1=peak_tp(p);
                        p2=peak_tp(p+1);
                        [min_v_cycle,min_tp_cycle]=min(est_val(p1:p2));
                        corrected_trough_values=[corrected_trough_values; min_v_cycle];
                    end
                    
                else
                    % less 2 peaks and no-nucleus
                    est_val=NaN*ones(length(tp),1);
                    est_b=NaN*ones(length(tp),1);
                end
            else
                % nucleus droplets
                est_val=NaN*ones(length(tp),1);
                est_b=NaN*ones(length(tp),1);
            end
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear6=est_val;']);
            tracks(n).Feat.est_linear_baseline=est_b;  
        end
        
        %% set average trough value to 1.0
        %disp(corrected_trough_values);
        mean_ctv=mean(corrected_trough_values);
        mean_estb=nanmean(corrected_baseline,2);
        %disp(mean_estb);
        for n =1:size(tracks)
            % get droplet and estimated data
            eval(['did=tracks(',num2str(n),').id;']);
            eval(['tmpSig=tracks(',num2str(n),').Feat.',colNames{1},';']);
            eval(['tp=tracks(',num2str(n),').Feat.t;']);
            %disp(n);
            %disp(size(mean_estb));
            tmpSig=tmpSig./mean_estb(tp(1):tp(end));
            eval(['tracks(',num2str(n),').Feat.',colNames{1},'_linear6=tmpSig;']);
        end
        
    end
end

%}                