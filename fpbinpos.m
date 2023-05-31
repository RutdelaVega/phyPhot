function [results] = fpbinpos(resultspath)

% Descriptives for position vs fiber photometry signal

% INPUTS:
%  resultspath: path where results was saved 

% OUTPUTS:
%  results: it saves all transformations of data (position-based) and
%  AUC/density calculations

% Figures are stored in same folder as results.mat

%% 1. Load results data:

load(resultspath, 'results')

% 2. Load data: wheel and FF

wheelnormcomplete = results.Behavior.Position.processed;
GCAMPcomplete = results.FP.Signals.raw.GCAMP;
DFFcomplete = results.FP.Signals.DFF;
DFFZcomplete = results.FP.Signals.DFFModZscore;
speedcomplete = results.Behavior.Speed;
Fs = results.FP.params.fs;

% Records with VR and no VR parts: 

% Wheel

if any(contains(fieldnames(results.Behavior), 'VR') == 1)

%wheelnorm = results.Behavior.Location.complete;

wheelnormVR = results.Behavior.PositionVR.processed;
wheelnormnoVR = results.Behavior.PositionnoVR.processed;

speedVR = results.Behavior.SpeedVR;
speednoVR = results.Behavior.SpeednoVR;

GCAMPVR = results.FP.SignalsVR.raw.GCAMP;
DFFVR = results.FP.SignalsVR.DFF;
DFFZVR = results.FP.SignalsVR.DFFModZscore;

GCAMPnoVR = results.FP.SignalsnoVR.raw.GCAMP;
DFFnoVR = results.FP.SignalsnoVR.DFF;
DFFZnoVR = results.FP.SignalsnoVR.DFFModZscore;

end


if any(contains(fieldnames(results.Behavior), 'VR') == 1)

wheeldata.wheelnormcomplete = wheelnormcomplete;
wheeldata.wheelnormVR = wheelnormVR;
wheeldata.wheelnormnoVR = wheelnormnoVR;

GCAMPdata.GCAMPcomplete = GCAMPcomplete;
GCAMPdata.GCAMPVR = GCAMPVR;
GCAMPdata.GCAMPnoVR = GCAMPnoVR;

DFFdata.DFFcomplete = DFFcomplete;
DFFdata.DFFVR = DFFVR;
DFFdata.DFFnoVR = DFFnoVR;

DFFZdata.DFFZcomplete = DFFZcomplete;
DFFZdata.DFFZVR = DFFZVR;
DFFZdata.DFFZnoVR = DFFZnoVR;

speeddata.speedcomplete = speedcomplete;
speeddata.speedVR = speedVR;
speeddata.speednoVR = speednoVR;

wheelnames = fieldnames(wheeldata);
GCAMPnames = fieldnames(GCAMPdata);
DFFnames = fieldnames(DFFdata);
DFFZnames = fieldnames(DFFZdata);
speednames = fieldnames(speeddata);

else

    wheeldata.wheelnormcomplete = wheelnormcomplete;
    GCAMPdata.GCAMPcomplete = GCAMPcomplete;
    DFFdata.DFFcomplete = DFFcomplete;
    DFFZdata.DFFZcomplete = DFFZcomplete;
    speeddata.speedcomplete = speedcomplete;

    wheelnames = fieldnames(wheeldata);
    GCAMPnames = fieldnames(GCAMPdata);
    DFFnames = fieldnames(DFFdata);
    DFFZnames = fieldnames(DFFZdata);
    speednames = fieldnames(speeddata);

end



% Other common variables: 

diameter = pi * 40;
speed_threshold = results.Behavior.params.MovThreshold; 


% Generate folders as a funtion of VR/noVR:

figpath = strcat(results.FP.path, '/Figures');
  for i = 1:length(wheelnames)
    mkdir(figpath, wheelnames{i})
    picpath{i} = strcat(figpath, '/' ,wheelnames{i});
  end




% Loop through all variables (complete, VR or no VR) dependent on if it's
% with VR or not. 

for i = 1:length(wheelnames)

wheelnorm = wheeldata.(wheelnames{i});
GCAMP = GCAMPdata.(GCAMPnames{i});
DFF = DFFdata.(DFFnames{i});
DFFZ = DFFZdata.(DFFZnames{i});
speed = speeddata.(speednames{i}); 

%% 2. Position binning


[pks, locs] = findpeaks(wheelnorm, 'MinPeakHeight', diameter - 2);

posbin = [0:(diameter*0.1):diameter];

locs = [1, locs]; locs = [locs, length(GCAMP)]; 


%% 3. Take FP values when the animal is moving (threshold saved in extractBeh):


%speed = speed(pre + 1 : end - post);

idx_thresh = find(speed >= speed_threshold);

% 3.1. Raw signal:

GCAMPval_thresh = GCAMP(speed >= speed_threshold);

% 3.2. DFF score:

DFFval_thresh = DFF(speed >= speed_threshold);

% 3.3. DFF Z-Score:

DFFZval_thresh = DFFZ(speed >= speed_threshold);

% 3.4. Wheel:

wheelval_thresh = wheelnorm(speed >= speed_threshold);

newsigvectorwheel = zeros(length(wheelnorm), 1).';
newsigvectorwheel(idx_thresh) = wheelval_thresh;
newsigvectorwheel(newsigvectorwheel == 0) = nan;

%% 4. Plot signals and position through time:

% 4.1. GCAMP (original raw signal):

fig1 = figure(1)
tiledlayout(3, 1)

nexttile
time = [0:length(GCAMP)-1]./Fs;
plot(time, speed, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
yline(speed_threshold, 'Color', 'red', 'LineWidth', 2)
ylabel('Speed (cm/s)')
hold off

nexttile

newsigvectorGCAMP = zeros(length(GCAMP), 1).';
newsigvectorGCAMP(idx_thresh) = GCAMPval_thresh;
newsigvectorGCAMP(newsigvectorGCAMP == 0) = nan;

plot(time, newsigvectorGCAMP, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on 
ylabel('GCAMP raw Signal')
hold off


nexttile

newposvectorGCAMP = zeros(length(GCAMP), 1).';
newposvectorGCAMP(idx_thresh) = wheelnorm(idx_thresh);
newposvectorGCAMP(newposvectorGCAMP == 0) = nan;

plot(time, newposvectorGCAMP, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
ylabel('Position (cm)')

xlabel('Time (s)')
hold off



% 4.2. DFF signal:

fig2 = figure(2)
tiledlayout(3, 1)

nexttile
time = [0:length(DFF)-1]./Fs;
plot(time, speed, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
yline(speed_threshold, 'Color', 'red', 'LineWidth', 2)
ylabel('Speed (cm/s)')
hold off

nexttile

newsigvectorDFF = zeros(length(DFF), 1).';
newsigvectorDFF(idx_thresh) = DFFval_thresh;
newsigvectorDFF(newsigvectorDFF == 0) = nan;

plot(time, newsigvectorDFF, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on 
ylabel('DFF Signal (%)')
hold off


nexttile

newposvectorDFF = zeros(length(DFF), 1).';
newposvectorDFF(idx_thresh) = wheelnorm(idx_thresh);
newposvectorDFF(newposvectorDFF == 0) = nan;

plot(time, newposvectorDFF, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
ylabel('Position (cm)')

xlabel('Time (s)')
hold off

% 4.3. DFFZ signal:

fig3 = figure(3)
tiledlayout(3, 1)

nexttile
time = [0:length(DFFZ)-1]./Fs;
plot(time, speed, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
yline(speed_threshold, 'Color', 'red', 'LineWidth', 2)
ylabel('Speed (cm/s)')
hold off

nexttile

newsigvectorDFFZ = zeros(length(DFFZ), 1).';
newsigvectorDFFZ(idx_thresh) = DFFZval_thresh;
newsigvectorDFFZ(newsigvectorDFFZ == 0) = nan;

plot(time, newsigvectorDFFZ, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on 
ylabel('/DeltaF/F (%)')
hold off


nexttile

newposvectorDFFZ = zeros(length(DFFZ), 1).';
newposvectorDFFZ(idx_thresh) = wheelnorm(idx_thresh);
newposvectorDFFZ(newposvectorDFFZ == 0) = nan;

plot(time, newposvectorDFFZ, 'LineWidth', 0.001)
xlim([min(time) max(time)])
hold on
ylabel('Position (cm)')

xlabel('Time (s)')
hold off


%% 5. Extract signal, position and AUC:

% Raw Signal, DFF Signal and DFFZ Signal:


for n = 1:(length(locs)-1)
        % Cut through each turn of the wheel and extract vectors
        temp = newsigvectorwheel(locs(n):locs(n+1)-1); 
        tempGCAMP = newsigvectorGCAMP(locs(n):locs(n+1)-1);
        tempDFF = newsigvectorDFF(locs(n):locs(n+1)-1);
        tempDFFZ = newsigvectorDFFZ(locs(n):locs(n+1)-1);

    for p = 1:(length(posbin)-1)
        % Within previous cut, identify and save values that belong to each
        % of the binned positions:

        posname = sprintf('Pos%i', p);
        sigwthreshGCAMP(n).(posname) = tempGCAMP((temp >= posbin(p)) & (temp < posbin(p+1)));
        sigwthreshDFF(n).(posname) = tempDFF((temp >= posbin(p)) & (temp < posbin(p+1)));
        sigwthreshDFFZ(n).(posname) = tempDFFZ((temp >= posbin(p)) & (temp < posbin(p+1)));

        bintimeGCAMP(n).(posname) = [p-1:1./(length(sigwthreshGCAMP(n).(posname))-1):p];
        bintimeDFF(n).(posname) = [p-1:1./(length(sigwthreshDFF(n).(posname))-1):p];
        bintimeDFFZ(n).(posname) = [p-1:1./(length(sigwthreshDFFZ(n).(posname))-1):p];

        % Calculate AUC:
        if length(bintimeGCAMP(n).(posname)) > 200

            AUCdataGCAMP(n, p) = trapz(bintimeGCAMP(n).(posname), sigwthreshGCAMP(n).(posname));
            AUCdataDFF(n, p) = trapz(bintimeDFF(n).(posname), sigwthreshDFF(n).(posname));
            AUCdataDFFZ(n, p) = trapz(bintimeDFFZ(n).(posname), sigwthreshDFFZ(n).(posname));
    
            densGCAMP(n, p) = AUCdataGCAMP(n, p)./(length(sigwthreshGCAMP(n).(posname))./Fs);
            densDFF(n, p) = AUCdataDFF(n, p)./(length(sigwthreshDFF(n).(posname))./Fs);
            densDFFZ(n, p) = AUCdataDFFZ(n, p)./(length(sigwthreshDFFZ(n).(posname))./Fs);

        else 
            AUCdataGCAMP(n, p) = nan;
            AUCdataDFF(n, p) = nan;
            AUCdataDFFZ(n, p) = nan;
    
            densGCAMP(n, p) = nan;
            densDFF(n, p) = nan;
            densDFFZ(n, p) = nan;
        end
    end

end


%% 6. Visualize signal through position:

fig4 = figure(4)

tiledlayout(3, 1)

nexttile

ylim([min(GCAMP) max(GCAMP)])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimeGCAMP(n).(posname), sigwthreshGCAMP(n).(posname), 'Color', 'black')
        hold on
        
    end
end
ylabel('GCAMP raw Signal')
hold off

nexttile

ylim([min(DFF) max(DFF)])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimeDFF(n).(posname), sigwthreshDFF(n).(posname), 'Color', 'black')
        hold on
        
    end
end
ylabel('DFF (%)')
hold off

nexttile

ylim([min(DFFZ) max(DFFZ)])
hold on
for n = 1:(length(locs)-1)
    for p = 1:(length(posbin)-1)
        posname = sprintf('Pos%i', p);
        plot(bintimeDFFZ(n).(posname), sigwthreshDFFZ(n).(posname), 'Color', 'black')
        hold on
        
    end
end
xlabel('Position')
ylabel('DFF Z-Score (%)')
hold off


%% 7. Plot AUC and density of AUC with different visualizations:

% Histogram:

AUCGCAMPmean = nanmean(AUCdataGCAMP).';
AUCDFFmean = nanmean(AUCdataDFF).';
AUCDFFZmean = nanmean(AUCdataDFFZ).';

posnames = fieldnames(bintimeGCAMP);
AUCdatatable = table(posnames, AUCGCAMPmean, AUCDFFmean, AUCDFFZmean);
AUCdatatable.posnames = categorical(AUCdatatable.posnames);

fig5 = figure(5)

tiledlayout(3, 1)

nexttile
bar(AUCdatatable.AUCGCAMPmean);
hold on 
ylim([max(AUCdatatable.AUCGCAMPmean - 3) max(AUCdatatable.AUCGCAMPmean)])
ylabel('AUC mean (GCAMP data)')
hold off

nexttile
bar(AUCdatatable.AUCDFFmean);
hold on
ylabel('AUC mean (DFF data)')
hold off

nexttile
bar(AUCdatatable.AUCDFFZmean);
hold on
ylabel('AUC mean (DFFZ data)')
hold off
 

fig13 = figure(13)

left_color = [0.92 0.69 0.12];right_color = [0 0 1];
set(fig13,'defaultAxesColorOrder',[left_color; right_color]);

tiledlayout(3, 1)

ix1 = nexttile

yyaxis(ix1, 'left')
bar(AUCdatatable.AUCGCAMPmean, "FaceColor",left_color);
hold on
ylim([0 max(AUCdatatable.AUCGCAMPmean)])
yyaxis(ix1, 'right')
bar(AUCdatatable.AUCGCAMPmean.', "FaceColor", right_color);
hold on 
ylim([0 max(AUCdatatable.AUCGCAMPmean)*2])
ylabel('AUC mean (GCAMP data)')

ix2 = nexttile
yyaxis(ix2, 'left')
bar(AUCdatatable.AUCDFFmean, "FaceColor",left_color);
hold on
ylim([0 max(AUCdatatable.AUCDFFmean)])
yyaxis(ix2, 'right')
bar(AUCdatatable.AUCDFFmean.', "FaceColor", right_color);
hold on 
ylim([0 max(AUCdatatable.AUCDFFmean)*2])
ylabel('AUC mean (DFF data)')


ix3 = nexttile
yyaxis(ix3, 'left')
bar(AUCdatatable.AUCDFFZmean,"FaceColor",left_color);
hold on
ylim([0 max(AUCdatatable.AUCDFFZmean)])
yyaxis(ix3, 'right')
bar(AUCdatatable.AUCDFFZmean.', "FaceColor",right_color);
hold on 
ylim([0 max(AUCdatatable.AUCDFFZmean)*2])
ylabel('AUC mean (DFFZ data)')



% Boxplots:

% AUC:

fig6 = figure(6)

tiledlayout(3, 1)


nexttile
boxchart(AUCdataGCAMP)
hold on
ylabel("AUC GCAMP")
title("AUC all Signals")
hold off

nexttile
boxchart(AUCdataDFF)
hold on
ylabel("AUC DFF")
hold off

nexttile
boxchart(AUCdataDFFZ)
hold on
ylabel("AUC DFFZ")
hold off


% Density:

%%% Boxplots


fig8 = figure(8)

tiledlayout(3, 1)

nexttile
boxchart(densGCAMP)
hold on
ylabel('GCAMP area density')
hold off

nexttile
boxchart(densDFF)
hold on
ylabel('DFF area density')
hold off


nexttile
boxchart(densDFFZ)
hold on
ylabel('DFFZ area density')
hold off

%%% Barplots:

% MEan density:
meandensGCAMP = nanmean(densGCAMP, 1);
meandensDFF = nanmean(densDFF, 1);
meandensDFFZ = nanmean(densDFFZ, 1);

fig9 = figure(9)

tiledlayout(3, 1)
set(fig9,'defaultAxesColorOrder',[left_color; right_color]);

ax1 = nexttile;
yyaxis(ax1, 'left')
bar(meandensGCAMP, "FaceColor",left_color)
hold on
ylabel('Mean Density (GCAMP)')
ylim([0 max(meandensGCAMP)])
yyaxis(ax1, 'right')
bar(meandensGCAMP,"FaceColor", right_color)
hold on
ylim([0 max(meandensGCAMP)*2])


ax2 = nexttile;
yyaxis(ax2, 'left')
bar(meandensDFF, "FaceColor",left_color)
hold on
ylabel('Mean Density (DFF)')
ylim([0 max(meandensDFF)])
yyaxis(ax2, 'right')
bar(meandensDFF, "FaceColor",right_color)
hold on
ylim([0 max(meandensDFF)*2])


ax3 = nexttile;
yyaxis(ax3, 'left')
bar(meandensDFFZ, "FaceColor",left_color)
hold on
ylabel('Mean Density (DFFZ)')
ylim([0 max(meandensDFFZ)])
yyaxis(ax3, 'right')
bar(meandensDFFZ, "FaceColor",right_color)
hold on
ylim([0 max(meandensDFFZ)*2])

% Median density:

mediandensGCAMP = nanmedian(densGCAMP, 1);
mediandensDFF = nanmedian(densDFF, 1);
mediandensDFFZ = nanmedian(densDFFZ, 1);

fig10 = figure(10)
set(fig10,'defaultAxesColorOrder',[left_color; right_color]);

tiledlayout(3, 1)

ax1 = nexttile;
yyaxis(ax1, 'left')
bar(mediandensGCAMP, "FaceColor",left_color)
hold on
ylabel('Median Density (GCAMP)')
ylim([0 max(mediandensGCAMP)])
yyaxis(ax1, 'right')
bar(mediandensGCAMP,"FaceColor", right_color)
hold on
ylim([0 max(mediandensGCAMP)*2])



ax2 = nexttile;
yyaxis(ax2, 'left')
bar(mediandensDFF, "FaceColor",left_color)
hold on
ylabel('Median Density (DFF)')
ylim([0 max(mediandensDFF)])
yyaxis(ax2, 'right')
bar(mediandensDFF,"FaceColor", right_color)
hold on
ylim([0 max(mediandensDFF)*2])


ax3 = nexttile;
yyaxis(ax3, 'left')
bar(mediandensDFFZ, "FaceColor",left_color)
hold on
ylabel('Mean Density (DFFZ)')
ylim([0 max(mediandensDFFZ)])
yyaxis(ax3, 'right')
bar(mediandensDFFZ, "FaceColor",right_color)
hold on
ylim([0 max(mediandensDFFZ)*2])



% Polar histogram:

fig7 = figure(7)

tiledlayout(1, 3)

nexttile
polarplot(AUCGCAMPmean)
title('Mean AUC GCAMP')

nexttile
polarplot(AUCDFFmean)
title('Mean AUC DFF')

nexttile
polarplot(AUCDFFZmean)
title('Mean AUC DFFZ')

%% 8. Save data in results.m:

% New vectors with speed threshold:

results.Analysis.wheelPosition.(wheelnames{i}).GCAMPvspeedthreshold = newsigvectorGCAMP;
results.Analysis.wheelPosition.(wheelnames{i}).DFFvspeedthreshold = newsigvectorDFF;
results.Analysis.wheelPosition.(wheelnames{i}).DFFZvspeedthreshold = newsigvectorDFFZ;


% Complete Data:

results.Analysis.wheelPosition.(wheelnames{i}).posvsGCAMP = sigwthreshGCAMP;
results.Analysis.wheelPosition.(wheelnames{i}).posvsDFF = sigwthreshDFF;
results.Analysis.wheelPosition.(wheelnames{i}).posvsDFFZ = sigwthreshDFF;
results.Analysis.wheelPosition.(wheelnames{i}).bintime = bintimeDFF;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.GCAMP = AUCdataGCAMP;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.GCAMPmean = AUCGCAMPmean;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.DFF = AUCdataDFF;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.DFFmean = AUCDFFmean;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.DFFZ = AUCdataDFFZ;
results.Analysis.wheelPosition.(wheelnames{i}).AUC.DFFZmean = AUCDFFZmean;
results.Analysis.wheelPosition.(wheelnames{i}).locs = locs; % Para la LFP de después 
results.Analysis.wheelPosition.(wheelnames{i}).wheel = wheelnormcomplete;

% Save Figures:


figpath1 = strcat(picpath{i}, '/GCAMPSpeedThreshVal.jpeg');
figpath2 = strcat(picpath{i}, '/DFFSpeedThreshVal.jpeg');
figpath3 = strcat(picpath{i}, '/DFFZSpeedThreshVal.jpeg');
figpath4 = strcat(picpath{i}, '/PosvsSig.jpeg');
figpath5 = strcat(picpath{i}, '/MeanAUChistogram.jpeg');
figpath6 = strcat(picpath{i}, '/AUCboxplot.jpeg');
figpath7 = strcat(picpath{i}, '/AUCpolar.jpeg');
figpath8 = strcat(picpath{i}, '/AreaDensity.jpeg');
figpath9 = strcat(picpath{i}, '/MeanDensity.jpeg');
figpath10 = strcat(picpath{i}, '/MedianDensity.jpeg');
figpath13 = strcat(picpath{i}, '/MeanAUChistogramscaled.jpeg');

saveas(fig1, figpath1);
saveas(fig2, figpath2);
saveas(fig3, figpath3);
saveas(fig4, figpath4);
saveas(fig5, figpath5);
saveas(fig6, figpath6);
saveas(fig7, figpath7);
saveas(fig8, figpath8);
saveas(fig9, figpath9);
saveas(fig10, figpath10);
saveas(fig13, figpath13);

resultspath = strcat(results.FP.path, '/Results.mat');
save(resultspath, 'results');


end % del bucle

end % de la función