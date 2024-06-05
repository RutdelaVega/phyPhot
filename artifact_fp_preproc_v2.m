
function [results] = artifact_fp_preproc_v2(TDTdata, varargin)

% Visualize your raw signals, correct artifacts if necessary and normalize
% data for further processing.

% INPUTS:
%   TDTdata (struct): fpdata.mat previously saved in fp2mat.m

% OPTIONAL INPUTS:
%   SampleBuffer (matrix): cuts signal at the beginning and end of recording.
%     Default is (5, 5) seconds. 
%   SmoothMethod (character): name of the method to observe tendencies in data more clearly. 
%     Default is 'Moving Average'. 
%   WindowSmooth (numeric): size of window for smoothing the signals
%     Default is 299. 
%   ZScoreMethod (cellstring): Z-transform of data for normalization. 
%     Default is 'Modified Median Z Score'.   

% OUTPUTS: 
%   Results: a struct containing all data (raw and processed) with all the
%   chosen parameters. 
% Data will be loaded to current workspace and saved in a results.mat file with all relevant information.
% Figures will be saved in the same folder where TDT data is stored and a 
% sub-folder "Figures" is created for that purpose. 


p = inputParser;
addParameter(p,'SampleBuffer', [5, 5],@ismatrix);
addParameter(p,'SmoothMethod', 'Moving Average', @ischar);
addParameter(p,'WindowSmooth', 299, @isnumeric);
addParameter(p,'ZScoreMethod', 'Modified Z-Score', @ischar);
addParameter(p, 'applyfilter', true, @islogical);
parse(p,varargin{:});
SampleBuffer = p.Results.SampleBuffer;
SmoothMethod = p.Results.SmoothMethod;
WindowSmooth = p.Results.WindowSmooth;
ZScoreMethod = p.Results.ZScoreMethod;
applyfilter = p.Results.applyfilter;


%% 1. Create a struct where all the data obtained with this function will be
% stored and create Figures Folder. Save figures path. 

results = struct();
mkdir(TDTdata.path, 'Figures')
figpath = strcat(TDTdata.path, '/Figures');


%% 2. Extract data

% 2.1. Extract single variables: 

Fs = TDTdata.streams.x465N.fs;

% Save to results

results.FP.params.fs = Fs;

% 2.2. Visualize raw signals:

GCAMPcheck = TDTdata.streams.x465N.data;
Isoscheck = TDTdata.streams.x405N.data;

Timecheck = (0:(length(GCAMPcheck)-1))./Fs;

fig1 = figure(1);
fig_width = 1800; % Specify the width of the figure window
fig_height = 400; % Specify the height of the figure window
fig_position = [50, 500, fig_width, fig_height]; % [left bottom width height]
set(gcf, 'Position', fig_position);

subplot(2, 1, 1) 
plot(Timecheck, GCAMPcheck, 'Color', [0,0.7,0.9]);
title('GCAMP raw signal')
% hold on;
subplot(2, 1, 2) 
plot(Timecheck, Isoscheck, 'Color', [0.4940 0.1840 0.5560]);
title('Isosbestic raw signal')
% hold off;
% legend('GCamP', 'Isosbestic');

fig1path = strcat(figpath, '/RawSignals.jpg');
saveas(fig1, fig1path); 

% 2.3. Optional input: Sample Buffer

% Save parameter to results

results.FP.params.SampleBuffer = SampleBuffer;

% Based on sample number (fs), how many samples do we need to cut:

pre = SampleBuffer(1);
post = SampleBuffer(2);
SBpre = round(pre .* Fs); 
SBpost = round(post .* Fs);

% Apply cut to data and create Time variable accordingly:

GCAMP = TDTdata.streams.x465N.data(SBpre + 1 : end - SBpost);
Isos = TDTdata.streams.x405N.data(SBpre + 1 : end - SBpost);

Time = (0:(length(GCAMP)-1))./Fs;


%% 3. Plot raw signal (to check for artifacts) and save result:


% Calculate signal quality: 

signalquality = corrcoef(GCAMP, Isos);
number = {'Corrcoef:' num2str(signalquality(1, 2))};

% Plot

fig2 = figure(2); 
fig_width = 1800; % Specify the width of the figure window
fig_height = 400; % Specify the height of the figure window
fig_position = [50, 500, fig_width, fig_height]; % [left bottom width height]
set(gcf, 'Position', fig_position);

subplot(2, 1, 1) 
plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
title('GCAMP raw signal')
% hold on;
subplot(2, 1, 2) 
plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
title('Isosbestic raw signal')
% hold off;
% legend('GCamP', 'Isosbestic');

fig2path = strcat(figpath, '/RawSignalswSB.jpg');
saveas(fig2, fig2path); 

results.FP.Signals.raw.GCAMP = GCAMP;
results.FP.Signals.raw.Isos = Isos;
results.FP.Signals.raw.Time = Time;



%% 3.0. Aplicar pass-band filter 15 Hz low-pass: 

% Filter
    
    if applyfilter == true
        fc=20;% cut off frequency
        fn=Fs; %nyquivst frequency = sample frequency/2;
        order = 2; %6th order filter, low pass
        [b14, a14]=butter(order,(fc/fn),'low');

        GCAMP = filtfilt(b14,a14,double(GCAMP)); 
        Isos = filtfilt(b14, a14, double(Isos));

        
        fig22 = figure(22); 
        set(gcf, 'Position', fig_position);
        subplot(2, 1, 1)
        plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
        title('GCAMP 15 Hz filtered')
        subplot(2, 1, 2)
        plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
        title('Isos 15 Hz filtered')

    end


%% 4. Bleaching correction 

% Julio Esparza's code: 

% Ask whether you want to correct for bleaching: 

correctbleaching = questdlg('Correction for bleaching?: ');

switch correctbleaching
    case 'Yes'
        correctbleaching = 1;
    case 'No'
        correctbleaching = 0;
end

if correctbleaching == 1
    [GCAMP, Isos, Time] = correctbleaching_FP(GCAMP, Isos, Time, Fs);
    results.FP.Signals.debleached.Time = Time;
    results.FP.Signals.debleached.GCAMP = GCAMP;
    results.FP.Signals.debleached.Isos = Isos;
end

%% 5. Prueba como Guppy para quedarnos con los trozos que nos interesan: 

keepsignal = questdlg('Would you like to remove part of the signal?');

switch keepsignal
    case 'Yes'
        keepsignal = 1;
    case 'No'
        keepsignal = 0;
end

signaltstokeep = [];

if keepsignal == 1

        fig31 = figure(31);
        set(gcf, 'Position', fig_position);
    subplot(2, 1, 1)
        plot(GCAMP, 'Color', [0,0.7,0.9]);
        title("Corrected GCAMP")
    subplot(2, 1, 2)
        plot(Isos, 'Color', [0.4940 0.1840 0.5560]);
        title("Corrected Isos")

    for ii = 1:100
       uiwait(msgbox('Select first onset, then offset of signal to keep', 'Instructions', "modal"));
       [x, ~, ~] = ginput(2);
       if x(1) < 0 
           x(1) = 1;
       elseif x(2) > length(GCAMP)
           x(2) = length(GCAMP);
       end
       signaltstokeep(ii, :) = [x(1) x(2)];

       keepmore = questdlg('Do you want to select more chunks to keep?:');

       switch keepmore
           case 'Yes'
               keepmore = 1;
           case 'No'
               keepmore = 0;
       end

       if keepmore == 1
           continue
       elseif keepmore == 0
           break
       end

    end
end

if keepsignal == 1

    if signaltstokeep(end, 2) > length(GCAMP)
        signaltstokeep(end, 2) = length(GCAMP);
    end

    signaltstokeep = int32(signaltstokeep);

    nanGCAMP = nan(1, length(GCAMP));
    nanIsos = nan(1, length(Isos));

    for ii = 1:size(signaltstokeep, 1)
        nanGCAMP(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = GCAMP(signaltstokeep(ii, 1):signaltstokeep(ii, 2));
        nanIsos(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = Isos(signaltstokeep(ii, 1):signaltstokeep(ii, 2));
    end


    fig32 = figure(32);
        set(gcf, 'Position', fig_position);
    subplot(2, 1, 1)
    plot(Time, nanGCAMP,'Color', [0,0.7,0.9])
        title("GCAMP to keep")
    subplot(2, 1, 2)
    plot(Time, nanIsos, 'Color', [0.4940 0.1840 0.5560])
        title("Isos to keep")

    GCAMP = nanGCAMP;
    Isos = nanIsos;
end

%% 5. Smooth signals:


    if SmoothMethod == 'Moving Average'
        method = 'moving';
        smwin = 299;
    elseif SmoothMethod == 'Lowess'
        method = 'lowess';
        smwin = 0.002;
    end
    
    if keepsignal == 1
        GCAMPsm = nan(1, length(GCAMP));
        Isossm = GCAMPsm;
        for ii = 1:size(signaltstokeep, 1)
            GCAMPsm(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = smooth(GCAMP(signaltstokeep(ii, 1):signaltstokeep(ii, 2)), smwin, method);
            Isossm(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = smooth(Isos(signaltstokeep(ii, 1):signaltstokeep(ii, 2)), smwin, method);
        end
        
    elseif keepsignal == 0
        GCAMPsm = smooth(GCAMP(~isnan(GCAMP)), smwin, method);
        Isossm = smooth(Isos(~isnan(Isos)), smwin, method);
    end

   fig5 = figure(5);
    set(gcf, 'Position', fig_position);
    subplot(2, 1, 1)
    plot(Time, GCAMPsm, 'Color', [0,0.7,0.9]);
    title("GCAMP smoothed")
    %     hold on;
    subplot(2, 1, 2)
    plot(Time, Isossm, 'Color', [0.4940 0.1840 0.5560]);
    %     hold off;
    %     legend('GCamP', 'Isosbestic');
    title("Isos smoothed")
    % Save data and plot:

    results.FP.params.filter.method = method;
    results.FP.params.filter.span = smwin;

    results.FP.Signals.smoothed.Time = Time;
    results.FP.Signals.smoothed.GCAMP = GCAMPsm.';
    results.FP.Signals.smoothed.Isos = Isossm.';

    fig5path = strcat(figpath, '/SmoothedSignals.jpg');
    saveas(fig5, fig5path);

%% 6. Normalize signals:

% 6.1.. DF/F of signal: 

if keepsignal == 1

DF = nan(1, length(GCAMP));

    for ii = 1:size(signaltstokeep, 1)
        
        bls = polyfit(Isossm(signaltstokeep(ii, 1):signaltstokeep(ii, 2)), GCAMPsm(signaltstokeep(ii, 1):signaltstokeep(ii, 2)), 1);
        yfit = bls(1).*Isossm(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) + bls(2);
        DeltaF = (GCAMPsm(signaltstokeep(ii, 1):signaltstokeep(ii, 2))-yfit) ./ (yfit);
        
        DF(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = DeltaF;
        
    end

fig6 = figure(6);
set(gcf, 'Position', fig_position);
% plot(Time, DeltaF, 'Color', [0,0.7,0.9])
plot( DeltaF, 'Color', [0,0.7,0.9])
title('Delta F of F')
    yline(0)
    DeltaF = DF;

    % Save data and figure DFF:

    results.FP.Signals.DFF = DeltaF.';

    fig6path = strcat(figpath, '/DFF.jpg');
    saveas(fig6, fig6path);


elseif keepsignal == 0

    bls = polyfit(Isossm(1:end), GCAMPsm(1:end), 1);
    yfit = bls(1).*Isossm + bls(2);
    DeltaF = (GCAMPsm(:)-yfit(:)) ./ (yfit(:));

    fig6 = figure(6);
    set(gcf, 'Position', fig_position);
    plot(Time, DeltaF, 'Color', [0,0.7,0.9])
    title('Delta F of F')

    % Save data and figure DFF:

    results.FP.Signals.DFF = DeltaF.';

    fig6path = strcat(figpath, '/DFF.jpg');
    saveas(fig6, fig6path);

end


% 6.2. Calculate Z-Score by the chosen method and plot the signal with
% selected methods and save data:


if keepsignal == 1
    if ZScoreMethod == "Z-Score"
        DFZscore = nan(1, length(GCAMP));
        for ii = 1:size(signaltstokeep, 1)
         DFZ = zscore(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2)));
         DFZscore(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = DFZ;
        end
        fig7 = figure(7);
        set(gcf, 'Position', fig_position);
        plot(Time, DFZscore, 'Color', [0,0.7,0.9]);
        hold on
        legend('Delta F Z-Score');
        hold off
        title('Delta F/F Z-Score')

        results.FP.Signals.DFFZscore = DFZscore.';

        fig7path = strcat(figpath, '/DFFZScore.jpg');
        saveas(fig7, fig7path);
    elseif ZScoreMethod == "Modified Z-Score"
        DFmodscore = nan(1, length(GCAMP));
        for ii = 1:size(signaltstokeep, 1)
           DFZmod = (0.6745.*(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2))-(median(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2))))))/mad(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2)), 1);
           DFmodscore(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = DFZmod;
        end
       fig7 = figure(7);
        set(gcf, 'Position', fig_position);
        plot(Time, DFmodscore, 'Color', [0,0.7,0.9]);
        hold on
        legend('Delta F Modified Z-Score');
        hold off
        title('Delta F/F Modified Z-Score')
        results.FP.Signals.DFFModZscore = DFmodscore.';

        fig7path = strcat(figpath, '/DFFModZScore.jpg');
        saveas(fig7, fig7path);
    end


elseif keepsignal == 0
    if ZScoreMethod == "Z-Score"
        
        DFZscore = zscore(DeltaF);
        fig7 = figure(7);
        set(gcf, 'Position', fig_position);
        plot(Time, DFZscore, 'Color', [0,0.7,0.9]);
        hold on
        legend('Delta F Z-Score');
        hold off
        title('Delta F/F Z-Score')


        results.FP.Signals.DFFZscore = DFZscore.';

        fig7path = strcat(figpath, '/DFFZScore.jpg');
        saveas(fig7, fig7path);


    elseif ZScoreMethod == "Modified Z-Score"

        DFmodscore = (0.6745.*(DeltaF-(median(DeltaF))))/mad(DeltaF, 1);
  fig7 = figure(7);
        set(gcf, 'Position', fig_position);
        plot(Time, DFmodscore, 'Color', [0,0.7,0.9]);
        hold on
        legend('Delta F Modified Z-Score');
        hold off
        title('Delta F/F Modified Z-Score')
        results.FP.Signals.DFFModZscore = DFmodscore.';

        fig7path = strcat(figpath, '/DFFModZScore.jpg');
        saveas(fig7, fig7path);
    end
end



%% 7. Save all parameters and data in a struct in chosen folder. 

resultspath = strcat(TDTdata.path, '/Results.mat');
results.FP.path = TDTdata.path;
save(resultspath, 'results');

end

