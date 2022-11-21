function plotStats(path,osciDat,dropdat,frameStep)
jitScale=0.4;
h=figure;
scatter(osciDat.periodPeak+rand(size(osciDat.periodPeak))*jitScale-jitScale/2,...
    osciDat.periodTrough+rand(size(osciDat.periodTrough))*jitScale-jitScale/2,[],...
    osciDat.cycleID,'.');
title('period comparison')
xlabel('period between peaks(frame)')
ylabel('period between troughs(frame)')
saveas(h,[path,'period_comparison.fig'])
close(h)

h=figure;
scatter(osciDat.timeLeft+rand(size(osciDat.periodPeak))*jitScale-jitScale/2,...
    osciDat.timeRight+rand(size(osciDat.periodTrough))*jitScale-jitScale/2,...
    [], osciDat.cycleID,'.');
title('period seperation, color label dropletID')
xlabel('peak to left trough(frame)')
ylabel('peak to right trough(frame)')
saveas(h,[path,'period_seperation.fig'])
close(h)

h=figure;
scatter(osciDat.ampLeft,...
    osciDat.ampRight,[], osciDat.cycleID,'.');
title('amp comparison, color label dropletID')
xlabel('peak to left trough amplitude')
ylabel('peak to right trough amplitude')
saveas(h,[path,'amp_comparison.fig'])
close(h)

h=figure;
scatter(osciDat.peakTime*frameStep+rand(size(osciDat.periodPeak))*jitScale*frameStep-jitScale*frameStep/2,...
    osciDat.periodPeak*frameStep+rand(size(osciDat.periodPeak))*jitScale*frameStep-jitScale*frameStep/2,...
    [],osciDat.dropID,'.')
title('time vs periodPeak, color label dropletID')
xlabel('Time')
ylabel('peak period')
saveas(h,[path,'time_peakPeriod.fig'])
close(h)

h=figure;
scatter(osciDat.peakTime*frameStep+rand(size(osciDat.periodPeak))*jitScale*frameStep-jitScale*frameStep/2,...
    osciDat.ampLeft,...
    [],osciDat.peakTime*frameStep,'.')
title('time vs ampLeft,color label peakTime')
xlabel('Time')
ylabel('amp left')
saveas(h,[path,'time_leftAmp.fig'])
close(h)

h=figure;
scatter(log10(osciDat.area),...
    osciDat.ampLeft,...
    [],osciDat.peakTime*frameStep,'.')
title('area vs ampLeft,,color label peakTime')
xlabel('log10(area),pixel^2')
ylabel('amp left')
saveas(h,[path,'area_leftAmp.fig'])
close(h)

h=figure;
scatter(log10(osciDat.area),...
    osciDat.periodPeak*frameStep+rand(size(osciDat.periodPeak))*jitScale*frameStep-jitScale*frameStep/2,...
    [],osciDat.peakTime*frameStep,'.')
title('area vs periodPeak')
xlabel('log10(area),pixel^2')
ylabel('periodPeak')
saveas(h,[path,'area_periodPeak.fig'])
close(h)
