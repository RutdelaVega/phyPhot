
function [results] = artifact_fp_preproc(TDTdata, varargin)

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
parse(p,varargin{:});
SampleBuffer = p.Results.SampleBuffer;
SmoothMethod = p.Results.SmoothMethod;
WindowSmooth = p.Results.WindowSmooth;
ZScoreMethod = p.Results.ZScoreMethod;


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
plot(Timecheck, GCAMPcheck, 'Color', [0,0.7,0.9]);
hold on;
plot(Timecheck, Isoscheck, 'Color', [0.4940 0.1840 0.5560]);
hold off;
legend('GCamP', 'Isosbestic');

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
plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
hold on
plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
hold off
legend('GCamP', 'Isosbestic');
text(Time(end), min(GCAMP), number,'Color','red','FontSize',8)

fig2path = strcat(figpath, '/RawSignalswSB.jpg');
saveas(fig2, fig2path); 

results.FP.Signals.raw.GCAMP = GCAMP;
results.FP.Signals.raw.Isos = Isos;
results.FP.Signals.raw.Time = Time;



% 3.0. Prueba con pass-band filter 15 Hz low-pass: 

   applyfilter = questdlg('Would you like to apply a 15 Hz low-pass filter on your data?');

        switch applyfilter
            case 'Yes'
               applyfilter = 1;
            case 'No'
               applyfilter = 0;
        end

% Filter
    
    if applyfilter == 1
        fc=20;% cut off frequency
        fn=Fs; %nyquivst frequency = sample frequency/2;
        order = 2; %6th order filter, low pass
        [b14, a14]=butter(order,(fc/fn),'low');

        GCAMP = filtfilt(b14,a14,double(GCAMP)); 
        Isos = filtfilt(b14, a14, double(Isos));

        
        fig22 = figure(22); 
plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
hold on;
plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
hold off;
legend('GCamP', 'Isosbestic');
title('signals with 15 Hz filter')

    end




% 3.1. Choose whether or not you want to correct for jump artifacts: 

correctartifact = questdlg('Would you like to correct for jump artifacts?');

switch correctartifact
    case 'Yes'
        correctartifact = 1;
    case 'No'
        correctartifact = 0;
end

if correctartifact == 1
    prompt = {'How many jump artifacts do you see? (Please count all artifacts separately, independent of signal):'};
    dlgtitle = 'Number of artifacts';
    dims = [1 35];
    definput = {'1'};
    art_quant = inputdlg(prompt, dlgtitle, dims, definput); % 

    art_quant = str2num(art_quant{1});

    results.FP.params.artifact.log = correctartifact;
    results.FP.params.artifact.quant = art_quant;

end


% 3.2. Select where to make the cut for jump artifacts (based on Julio
% Esparza's function: Fiber_photometry_GCaMPf.m):

if correctartifact == 1

      for nn = 1:art_quant

        
        signaltocorrect = questdlg('Which signal do you want to correct?', ...
        	'Select Signal', ...
        	'GCAMP','Isos', 'Isos');
        % Handle response
        switch signaltocorrect
            case 'GCAMP'
                signaltocorrect = 1;
            case 'Isos'
                signaltocorrect = 2;
        end
        
        if signaltocorrect == 1
            uiwait(msgbox('Select first left side, then right side of jump', 'Instructions', "modal"));
            [x, y, button] = ginput(2);

            leftidx = floor((x(1)).*Fs);
            rightidx = round((x(2)).*Fs);
        
            GCAMPtocorrect = GCAMP(leftidx:rightidx);
            Timetocorrect = Time(leftidx:rightidx);
        elseif signaltocorrect == 2
            uiwait(msgbox('Select first left side, then right side of jump', 'Instructions', "modal"));
            [x, y, button] = ginput(2);

            leftidx = floor((x(1)).*Fs);
            rightidx = round((x(2)).*Fs);
        
            Isostocorrect = Isos(leftidx:rightidx);
            Timetocorrect = Time(leftidx:rightidx);
        end



    % For better precision, it shows the selected jump again amplified

        if signaltocorrect == 1
            plot(Timetocorrect, GCAMPtocorrect);
            uiwait(msgbox('Select first left side, then right side of jump', 'Instructions', "modal"));
            [x, y, button] = ginput(2);

            leftidx = floor((x(1)).*Fs);
            rightidx = round((x(2)).*Fs);

            GCAMPtocorrect = GCAMP(leftidx:rightidx); 
            Timetocorrect = Time(leftidx:rightidx);
        elseif signaltocorrect == 2   
            plot(Timetocorrect, Isostocorrect)
            uiwait(msgbox('Select first left side, then right side of jump', 'Instructions', "modal"));
            [x, ~, ~] = ginput(2);

            leftidx = floor((x(1)).*Fs);
            rightidx = round((x(2)).*Fs);

            Isostocorrect = GCAMP(leftidx:rightidx); 
            Timetocorrect = Time(leftidx:rightidx);
        end

    if correctartifact == 1
        artifactnumber = sprintf('art%i', nn);
        if signaltocorrect == 1
        results.FP.params.artifact(nn).(artifactnumber) = [Timetocorrect(:), GCAMPtocorrect(:)];
        elseif signaltocorrect == 2
        results.FP.params.artifact(nn).(artifactnumber) = [Timetocorrect(:), Isostocorrect(:)];
        end
    end

% 3.3. Interpolate jump values and plot again with corrected signal:

        
        if signaltocorrect == 1
            GCAMP(leftidx+1:rightidx-1) = 0;
            yy = interp1(Time(1:end), GCAMP(1:end), GCAMP(leftidx+1:rightidx-1), 'spline');
            GCAMP(leftidx+1:rightidx-1) = yy;
            
            % Standardize values based on median:
            
            mRef = median(GCAMP(1:leftidx));
            stdRef = std(GCAMP(1:leftidx));
            mSig = median(GCAMP(rightidx:end));
            stdSig = std(GCAMP(rightidx:end));
            newGCAMP = ((GCAMP(rightidx:end) - mSig)/stdSig) .* stdRef + mRef;
        
            GCAMP = [GCAMP(1:rightidx-1), newGCAMP];
            
            % Save times and signal that needs to be scalated after artifacts (partial signal after artifact):
                       
%             results.FP.params.artifact.scalated(nn).Time = Time(rightidx:end);
%             results.FP.params.artifact.scalated(nn).signal = newGCAMP;
            
        elseif signaltocorrect == 2

            Isos(leftidx+1:rightidx-1) = 0;
            yy = interp1(Time(1:end), Isos(1:end), Isos(leftidx+1:rightidx-1), 'spline');
            Isos(leftidx+1:rightidx-1) = yy;
            
            % Standardize values based on median:
            
            mRef = median(Isos(1:leftidx));
            stdRef = std(Isos(1:leftidx));
            mSig = median(Isos(rightidx:end));
            stdSig = std(Isos(rightidx:end));
            newIsos = ((Isos(rightidx:end) - mSig)/stdSig) .* stdRef + mRef;
        
            Isos = [Isos(1:rightidx-1), newIsos];
            
            % Save times and signal that needs to be scalated after artifacts (partial signal after artifact):
            
            
%             results.FP.params.artifact.scalated(nn).Time = Time(rightidx:end);
%             results.FP.params.artifact.scalated(nn).signal = newIsos;

        end

                % Plot signal again with corrected signals:

        fig3 = figure(3);
        plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
        hold on
        plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
        hold off
        legend('GCamP', 'Isosbestic');
        title('Corrected Signals')


        % Save figure 3:

        fig3path = strcat(figpath, '/CorrectedSignals.jpg');
        saveas(fig3, fig3path);


      end

end

% Save corrected data to results (complete signal) 

if correctartifact == 1
    results.FP.Signals.corrected.Time = Time;
    results.FP.Signals.corrected.GCAMP = GCAMP;
    results.FP.Signals.corrected.Isos = Isos;
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
    GCAMP = GCAMP.';
    Isos = Isos.';
    Time = Time.';

    div_len = 40; 
    L = size(GCAMP,1); 
    div_len = floor(div_len*Fs);
    num_div = floor((L-floor(Fs*10))/div_len); 
    last_div = (L-floor(Fs*10)) - num_div*div_len;
    [~, ind] = min(GCAMP(1:floor(10*Fs)));
    [~, indIso] = min(Isos(1:floor(10*Fs)));

    temp_t = Time(ind);
    temp_s = GCAMP(ind);
    temp_i = Isos(indIso);
    tempGCAMP = GCAMP(floor(Fs*10)+1:end);
    tempIsos = Isos(floor(Fs*10)+1:end);
    
    adj = floor(10*Fs);

    for ii = 0:num_div-1
    a = tempGCAMP((ii*div_len)+1:(ii+1)*div_len);
    b = tempIsos((ii*div_len)+1:(ii+1)*div_len);
    a = sort(a);
    b = sort(b);
    a = a(1:round(Fs));
    b = b(1:round(Fs));
    for jj = 1:round(Fs)
        temp2GCAMP = find(tempGCAMP==a(jj));
        temp2Isos = find(tempIsos == b(jj));

        if ~isempty(temp2GCAMP)
            tempGCAMP(temp2GCAMP(1)) = +1000;
            ind(jj,1) = temp2GCAMP(1)+adj;
        else
            ind(jj,1) = NaN;
        end

        if ~isempty(temp2Isos)
            tempIsos(temp2Isos(1)) = + 1000;
            indIso(jj, 1) = temp2Isos(1)+adj;
        else
            indIso(jj, 1) = NaN;
        end

    end
    ind(isnan(ind)) = [];
    indIso(isnan(indIso)) = [];
    temp_t = [temp_t; Time(ind)];
    temp_s = [temp_s; GCAMP(ind)];
    temp_i = [temp_i; Isos(indIso)];
    ind = [];
    indIso = [];
    end 

    if last_div>=floor(25*Fs)
    a = tempGCAMP((ii+1)*div_len+1:end);
    b = tempIsos((ii+1)*div_len+1:end);
    a = sort(a);
    b = sort(b);
    a = a(1:round(Fs));
    b = b(1:round(Fs));
    for jj = 1:round(Fs)
        temp2GCAMP = find(tempGCAMP==a(jj));
        temp2Isos = find(tempIsos == b(jj));
        if ~isempty(temp2GCAMP)
            tempGCAMP(temp2GCAMP(1)) = + 1000;
            ind(jj,1) = temp2GCAMP(1)+adj;
        else
            ind(jj,1) = NaN;
        end
        if ~isempty(temp2Isos)
            tempIsos(temp2Isos(1)) = + 1000;
            indIso(jj, 1) = temp2Isos(1)+adj;
        else
            indIso(jj, 1) = NaN;
        end
    end
    ind(isnan(ind)) = [];
    temp_t = [temp_t; Time(ind)];
    temp_s = [temp_s; GCAMP(ind)];
    ind = [];
    end 

% GCAMP:

yy = interp1(temp_t,temp_s,Time, 'linear');
%fix NaNs at the beginning and at the end
a = isnan(yy);
ind2 = find(a==0);
st = ind2(1)-1;
fn = ind2(end)+1;
yy(1:st) = yy(st+1);
yy(fn:end) = yy(fn-1);
yy = movmean(yy, div_len);
yy2 = yy;

options = optimset('MaxFunEvals',10000);
fun2 = @(x) sum((yy2 - (x(1)*Time + x(2))).^2);
x0 = [0 min(GCAMP)];
bestx = fminsearch(fun2,double(x0), options); 
yy3 = bestx(1)*Time + bestx(2);

yy(:,1) = max([real(yy2), real(yy3)],[],2);
yy(:,2) = min([real(yy2), real(yy3)],[],2); 
yy = (0.75.*yy(:,1) + 0.25.*yy(:,2));
basal = real(yy);

% Isos:

yyI = interp1(temp_t,temp_i,Time, 'linear');
%fix NaNs at the beginning and at the end
aI = isnan(yyI);
ind2I = find(aI==0);
stI = ind2I(1)-1;
fnI = ind2I(end)+1;
yyI(1:stI) = yyI(stI+1);
yyI(fnI:end) = yyI(fnI-1);
yyI = movmean(yyI, div_len);
yy2I = yyI;

options = optimset('MaxFunEvals',10000);
fun2I = @(x) sum((yy2I - (x(1)*Time + x(2))).^2);
x0I = [0 min(Isos)];
bestxI = fminsearch(fun2I,double(x0I), options); 
yy3I = bestxI(1)*Time + bestxI(2);

yyI(:,1) = max([real(yy2I), real(yy3I)],[],2);
yyI(:,2) = min([real(yy2I), real(yy3I)],[],2); 
yyI = (0.75.*yyI(:,1) + 0.25.*yyI(:,2));
basalI = real(yyI);

GCAMP = GCAMP - basal + min(basal);
Isos = Isos - basalI + min(basalI);

fig4 = figure(4); 
    plot(Time, GCAMP, 'Color', [0,0.7,0.9]);
    hold on;
    plot(Time, Isos, 'Color', [0.4940 0.1840 0.5560]);
    hold off;
    legend('GCamP', 'Isosbestic');
    title('Debleached Signals')

    GCAMP = GCAMP.';
    Isos = Isos.';
    Time = Time.';
end

if correctbleaching == 1

    results.FP.Signals.debleached.Time = Time;
    results.FP.Signals.debleached.GCAMP = GCAMP;
    results.FP.Signals.debleached.Isos = Isos;

end

%% Prueba como Guppy para quedarnos con los trozos que nos interesan: 

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
        plot(GCAMP, 'Color', [0,0.7,0.9]);
        hold on
        plot(Isos, 'Color', [0.4940 0.1840 0.5560]);
        hold off
        legend('GCamP', 'Isosbestic');
        title('Corrected Signals')



    for ii = 1:100
       uiwait(msgbox('Select first onset, then offset of signal to keep', 'Instructions', "modal"));
       [x, ~, ~] = ginput(2);
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
    plot(Time, nanGCAMP,'Color', [0,0.7,0.9])
    hold on
    plot(Time, nanIsos, 'Color', [0.4940 0.1840 0.5560])
    hold off
    legend('GCamP', 'Isosbestic');
    title('Signals to keep')

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
    plot(Time, GCAMPsm, 'Color', [0,0.7,0.9]);
    hold on;
    plot(Time, Isossm, 'Color', [0.4940 0.1840 0.5560]);
    hold off;
    legend('GCamP', 'Isosbestic');
    title('Smoothed Signals')

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
    plot(Time, DF, 'Color', [0,0.7,0.9])
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
           DFZmod = (0.6745.*(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2))-(median(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2))))))/mad(DeltaF(signaltstokeep(ii, 1):signaltstokeep(ii, 2)));
           DFmodscore(signaltstokeep(ii, 1):signaltstokeep(ii, 2)) = DFZmod;
        end
        fig7 = figure(7);
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
        plot(Time, DFZscore, 'Color', [0,0.7,0.9]);
        hold on
        legend('Delta F Z-Score');
        hold off
        title('Delta F/F Z-Score')

        results.FP.Signals.DFFZscore = DFZscore.';

        fig7path = strcat(figpath, '/DFFZScore.jpg');
        saveas(fig7, fig7path);


    elseif ZScoreMethod == "Modified Z-Score"

        DFmodscore = (0.6745.*(DeltaF-(median(DeltaF))))/mad(DeltaF);
        fig7 = figure(7);
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

