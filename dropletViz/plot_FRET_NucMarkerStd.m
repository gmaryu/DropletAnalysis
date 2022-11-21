% load dataset
%load('Pos12tracking_20191116T022743.mat')
tracks=sumTrackTable{3};
%% add feild for orifinal id info
for i = 1:size(tracks, 1)
    %disp(i)
    tracks(i).id=i;
end
%% 
% thresholding by tracking frame number (75%?)
fnum = 300;
thres = 0.9;

for i =1:size(tracks, 1)
    if tracks(i).Len > fnum * thres
        % FRET ratio and YFP intensity
        %FRETratio = tracks(i).Feat.FRET ./ tracks(i).Feat.CFP;
        FRETratio_n = tracks(i).Feat.FRET_nuc ./ tracks(i).Feat.CFP_nuc;
        FRETratio_c = tracks(i).Feat.FRET_cyto ./ tracks(i).Feat.CFP_cyto;
        FRETratio_w = tracks(i).Feat.FRET_whole ./ tracks(i).Feat.CFP_whole;
        
        FYratio_n = tracks(i).Feat.FRET_nuc ./ tracks(i).Feat.YFP_nuc;
        FYratio_c = tracks(i).Feat.FRET_cyto ./ tracks(i).Feat.YFP_cyto;
        FYratio_w = tracks(i).Feat.FRET_whole ./ tracks(i).Feat.YFP_whole;
        
        % detrending ratio by YFP
        %nrm_FRETratio = FRETratio ./ YFP;
        %nrm_scl_FRETratio = nrm_FRETratio/min(nrm_FRETratio); %baseline rises to 1 

        %{
        % detrending ratio by ratio for peakdetection
        [p,s,mu] = polyfit((1:numel(FRETratio))',FRETratio,6);
        f_y = polyval(p,(1:numel(FRETratio))',[],mu);
        FRETratio_detrend = FRETratio - f_y; 
        %}
        
        % plot figure
        %figure('Position', [100, 100, 1024, 1200]);
        figure();
        hold on
        title(['pGM-020 #',num2str(tracks(i).id)]);
        xlabel('Time (frame)');
        
        % plot FRETratio
        yyaxis left;
        ylabel('FRET/CFP ratio');
        plot(FRETratio_n, 'g');
        plot(FRETratio_c, 'b');
        
        % plot NLS-mCherry_stdev
        yyaxis right;
        ylabel('NLS-mCherry-NLS Std');
        plot(tracks(i).Feat.mCherrystdInt, 'r');
        
        legend('FRET-nuc','FRET-cyto','MarkerStd');
        % save fig as png
        fname = ['TrackID_',num2str(tracks(i).id),'.png'];
        disp(fname);
        saveas(gcf, fname);
        close;
    end
end