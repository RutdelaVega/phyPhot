
function [results] = PETHZscore(resultspath, method, Pre, Post, bin, varargin)

% Based on TrialPETHZscore.m function from pMat (BarkerLab; Barker et al, 2021) 

% INPUTS:
%  resultspath: path where file results from previous function was saved
%  method: method for creating heatmap (take DFF directly or calculate
%  based on baseline)
%  pre: time before event to visualize
%  post: time after start of event to visualize
%  bin: sampling size per second

% OPTIONAL INPUTS: (in the case of baseline method)
%  windowsamplestart: time where baseline starts
%  windowduration: baseline window duration

% OUTPUTS:
%  results: struct with all the data


p = inputParser;
addParameter(p,'WinSampleStart', -3,@isnumeric);
addParameter(p,'WinDuration', 3, @isnumeric);
parse(p,varargin{:});
WinSampleStart = p.Results.WinSampleStart;
WinDuration = p.Results.WinDuration;


%% 1. Load results.mat:

load(resultspath, 'results')




%% 2. Check if "corrected" fieldname exists and extract variables: 

Fs = results.FP.params.fs;
BehFs = results.Behavior.Fs;
EventTS = results.Behavior.Event.Time./BehFs; 

iscorrected = isfield(results.FP.Signals, 'corrected');

if iscorrected == 1
    Ch490 = results.FP.Signals.corrected.GCAMP;
    Ch405 = results.FP.Signals.corrected.Isos;
    Ts = ([0:(length(Ch490)-1)])./Fs;
elseif iscorrected == 0
    Ch490 = results.FP.Signals.raw.GCAMP;
    Ch405 = results.FP.Signals.raw.Isos;
    Ts = ([0:(length(Ch490)-1)])./Fs;
end

DFFZ = double(results.FP.Signals.DFFModZscore);
DFF = results.FP.Signals.DFF;


%% 3. Apply filter if necessary:

  figure(1)
    plot(Ts, DFFZ);
    grid on;
    hold on;
    xlabel('Time');
    ylabel('\DeltaF/F (%)');
    hold off


   applyfilter = questdlg('Would you like to apply a 0.01 Hz high-pass filter on your data?');

        switch applyfilter
            case 'Yes'
               applyfilter = 1;
            case 'No'
               applyfilter = 0;
        end

% Filter
    
    if applyfilter == 1
        fc=0.01;% cut off frequency
        fn=Fs; %nyquivst frequency = sample frequency/2;
        order = 2; %6th order filter, high pass
        [b14, a14]=butter(order,(fc/fn),'high');

        DFFZ = filtfilt(b14,a14,DFFZ); 
    end





%% 4. Use parameters to calculate PETH, AUC and Peaks (same as in pMat):

 bin = round(bin*Fs);
 PreWind = round(Pre*Fs);
 PostWind = round(Post*Fs);
            
log = (EventTS(:, 2)-EventTS(:, 1)) >= Post;
eventstokeep = sum(log);
totalevents = length(EventTS);

info = sprintf('There are %i events (out of %i) that are at least as long as the chosen "Event Post Time" and that will be used for further analysis. ', eventstokeep, totalevents);
disp(info)

 
%% 5. Eliminate events with overlapping windows.

    tmp=[];CurrEvent=EventTS(1);
    %REMOVE -1 HERE
    for i=1:length(EventTS)-1
        if EventTS(i+1)-CurrEvent>Pre+Post
            tmp(end+1,1)=CurrEvent;
            CurrEvent=EventTS(i+1);
        else
        end
    end

    tmp(end+1,1)=EventTS(length(EventTS),1);
    EventTS=tmp;

%% 6. Find Time=0 for the event within the photometry data

    CSidx=[];
    for i=1:length(EventTS)
        [MinVal, CSidx(i,1)]=min(abs(Ts(:,:)-EventTS(i))); 
    end 

%% 7. Eliminate events too close to the beginning:

for i = 1:length(CSidx)
    if (CSidx(i)-PreWind) <= 0 
        CSidx(i) = nan;
    elseif (CSidx(i) + PostWind) > length(Ts)
        CSidx(i) = nan;
    end
end

CSidx = CSidx(~isnan(CSidx));

%% 8. Select method for computing and creating heatmaps:

if method == 'DFF'
    method = 1;
elseif method == 'baseline'
    method = 2;
end

if method == 1
    for i = 1:length(CSidx)
        CSTS = (-PreWind:PostWind)./Fs;
        DF_ZScore(:, i) = DFFZ((CSidx(i)-PreWind):(CSidx(i)+PostWind));
        DF_F(:, i) = DFF((CSidx(i)-PreWind):(CSidx(i)+PostWind));
    end
end


if method == 2

% 8.1. Baseline Sampling Window Params:

BaselineWind = WinSampleStart; BaselineWind2 = WinDuration;

BaselineWind=round(BaselineWind*Fs); 
BaselineWind2=round(BaselineWind2.*Fs); 

% 8.2. Obtain the DeltaF/F for each event window

    CS405=[];CSTS=[];CS490=[] ;
    Ch490=Ch490';
    Ch405=Ch405';
    counter=1;
    for i=1:length(CSidx)
        if CSidx(i)-(BaselineWind+BaselineWind2)<=0 || CSidx(i)+PostWind > length(Ts)
        else
            CSTS=(-PreWind:PostWind)./Fs;
            CSBL(1,:)=Ch490((CSidx(i)-(BaselineWind+BaselineWind2)):(CSidx(i)-(BaselineWind))); 
            CSBL2(1,:)=Ch405((CSidx(i)-(BaselineWind+BaselineWind2)):(CSidx(i)-(BaselineWind)));
            if i>length(CSidx)
             break
            elseif i<=length(CSidx) 

            CS405(1,:)=Ch405((CSidx(i)-PreWind):(CSidx(i)+PostWind)); 
            CS490(1,:)=Ch490((CSidx(i)-PreWind):(CSidx(i)+PostWind)); 
            end

            %Smooth to eliminate high frequency noise.
            F490=smooth(CS490,0.002,'lowess');  %Was 0.002- DJB
            F405=smooth(CS405,0.002,'lowess'); 
            F490CSBL=smooth(CSBL,0.002,'lowess');  %Was 0.002- DJB
            F405CSBL=smooth(CSBL2,0.002,'lowess');


            
            %Scale and fit data
            bls=polyfit(F405(1:end),F490(1:end),1);
            blsbase=polyfit(F405CSBL(1:end),F490CSBL(1:end),1);
            Y_Fit=bls(1).*F405+bls(2);
            Y_Fit_base=blsbase(1).*F405CSBL+blsbase(2);

            %Center data and generate Delta F/F (DF_F) by dividing
            %event window by median of baseline window.
            DF_Event(:,i)=(F490(:,1)-Y_Fit(:,1)); 
            DF_F(:,i)=DF_Event(:,i)./(Y_Fit); 
            DF_F(:,i)=DF_F(:,i)-DF_F(1,i);
            
            DF_Event(:,i)=(F490(:,1)-Y_Fit(:,1));
            DF_Event(:,i)=DF_Event(:,i)-DF_Event(1,i); 
            DF_Base(:,i)=(F490CSBL(:,1)-Y_Fit_base(:,1));
            DF_Base(:,i)=DF_Base(:,i)-DF_Base(1,i); 
            DF_ZScore(:,counter)=(DF_Event(:,i)-median(DF_Base(:,i)))./mad(DF_Base(:,i)); 
            counter=counter+1;
            %%% Clearing variables to reset for the next trial        
            clear CS405 CS490 F490 F405 bls Y_Fit

        end 
    end 
end 

%% 9. Plot PETH:

% DFF:
DF_F=DF_F';
tmp=[];tmp2=[];
    for i=1:bin:length(CSTS)
        if i+bin>length(CSTS)
        tmp(1,end+1)=median(CSTS(i:end));
        tmp2(:,end+1)=median(DF_F(:,i:end),2);
        else
        tmp(1,end+1)=median(CSTS(i:i+bin));
        tmp2(:,end+1)=median(DF_F(:,i:i+bin),2);
        end
    end
    CSTS3=tmp;
    DF_F=tmp2;

CSTrace2=(mean(DF_F,1)*100);


% ZScore DFF:
DF_ZScore=DF_ZScore';

tmp=[];tmp2=[];
    for i=1:bin:length(CSTS)
        if i+bin>length(CSTS)
        tmp(1,end+1)=median(CSTS(i:end));
        tmp2(:,end+1)=median(DF_ZScore(:,i:end),2);
        else
        tmp(1,end+1)=median(CSTS(i:i+bin));
        tmp2(:,end+1)=median(DF_ZScore(:,i:i+bin),2);
        end
    end
    CSTS2=tmp;
    DF_ZScore=tmp2;


fig1 = figure(1);
subplot(2,1,1)
imagesc(CSTS2,[1:size(DF_ZScore,1)],(DF_ZScore)); 
ylabel('Trial #','fontsize', 18)
title('Trial Based Z-Score PETH','fontsize', 18)

xlim([-Pre Post])
set(gca,'fontsize',18);
CSTrace1=(mean(DF_ZScore,1)); 
CSmax=max(max(CSTrace1)); 
CSmin=min(min(CSTrace1));

if CSmax<0
    CSmax=0.1;
end

TrialData=DF_ZScore;


subplot(2,1,2), 
plot(CSTS2,CSTrace1,'LineWidth',3) 
title('Z-Score PETH', 'Fontsize', 18);
xlim([-Pre Post])
xlabel('Time (s)', 'FontSize',18)
ylabel('Z-score','fontsize', 18);
ylim([CSmin CSmax*1.25]);
hold on
plot([0,0],[CSmin,CSmax*1.25],':r')
hold off


%% 10. Save data to calculate AUC and peaks

CSTrace1=CSTrace1'; CSTrace2=CSTrace2';

FinalData(:,1)=CSTS2';FinalData(:,2)=CSTrace1; FinalData(:,3)=CSTrace2; 

%% 11. Compute AUC and Peaks: 

AUC=[];MaxValue=[];

prompt3 = {'First Interval Start(s):', 'Second Interval Start (s):', 'Third Interval Start (s):'};

dlgtitle3 = 'Set AUC/Peak Time Window';
dims3 = [1 60];
definput3 = {'-2', '0', '2'}; 

SB3 = inputdlg(prompt3, dlgtitle3, dims3, definput3); 

first = str2double(SB3{1});
second = str2double(SB3{2});
third = str2double(SB3{3});
fourth = third + second-first; % En principio sí está bien así, lo digo por los signos menos y buscar la diferencia.. es para encontrar el intervalo que se ha usado y automáticamente calcular hasta donde llega el último intervalo

while (first < -Pre || fourth > Post) || ((third-second) ~= (second-first))

    mes = sprintf('First Interval: %i. Second Interval: %i. Third Interval: %i', -Pre, -Pre+1, -Pre+2);
    mes2 = 'Values must be within chosen event window or intervals need to be the same size. Please select new parameters (e.g.):';
    warnmes = strcat(mes2, mes);
    waitfor(warndlg(warnmes, 'Select New AUC/Peaks Parameters'))

    prompt3 = {'First Interval Start(s):', 'Second Interval Start (s):', 'Third Interval Start (s):'};

    dlgtitle3 = 'Set AUC/Peak Time Window';
    dims3 = [1 60];
    definput3 = {'-2', '0', '2'}; 

    SB3 = inputdlg(prompt3, dlgtitle3, dims3, definput3); 

    first = str2double(SB3{1});
    second = str2double(SB3{2});
    third = str2double(SB3{3});
    fourth = third + second-first;

end



% Take AUC Between specific timestamps from PETH data.

AUC(1,1)=trapz(FinalData(FinalData(:,1)>=first & FinalData(:,1)<=second,2)); % Está comprobando los índices de la segunda columna que estén por debajo de ese momento en el tiempo y extrae los valores de la DF_ZScore 
AUC(1,2)=trapz(FinalData(FinalData(:,1)>=second & FinalData(:,1)<=third,2));
AUC(1,3)=trapz(FinalData(FinalData(:,1)>=third & FinalData(:,1)<=fourth,2));
% 

% Take peak between specific timestamps from PETH data.


MaxValue(1,1)=max(FinalData(FinalData(:,1)>=first & FinalData(:,1)<=second,2)); % Está comprobando los índices de la segunda columna que estén por debajo de ese momento en el tiempo y extrae los valores de la DF_ZScore 
MaxValue(1,2)=max(FinalData(FinalData(:,1)>=second & FinalData(:,1)<=third,2));
MaxValue(1,3)=max(FinalData(FinalData(:,1)>=third & FinalData(:,1)<=fourth,2));


%% 12. Save data for later stats:   

results.Analysis.PETH.ZScore = DF_ZScore; 
results.Analysis.PETH.ZScoremean = CSTrace1.'; 
results.Analysis.PETH.DFF = DF_F; 
results.Analysis.PETH.DFFmean = CSTrace2.';

results.Analysis.AUC = AUC;
results.Analysis.Peaks = MaxValue;

results.Analysis.Params.PreWind = Pre;
results.Analysis.Params.PostWind = Post;
results.Analysis.Params.BinConstant = bin;

if method == 2
    results.Analysis.Params.SampleStartWind = BaselineWind;
    results.Analysis.Params.DurWind = BaselineWind2;
end

% Create Data Folder:

mkdir(results.FP.path, 'Data')
figpath = strcat(results.FP.path, '/Data');

% Save Data in same file: 

resultspath = strcat(results.FP.path, '/Results.mat');
save(resultspath, 'results');

% Save as CSV or XLSX:

colnames1 = {'Event Time Window', 'Mean DF_ZScore for Event', 'Mean DF_F for Event' };
FinalData = array2table(FinalData, 'VariableNames', colnames1);
xlsx1 = strcat(results.FP.path, '/Data/EventDF.xlsx');

colnames2 = {'First Interval', 'Second Interval', 'Third Interval'};
AUC = array2table(AUC, 'VariableNames', colnames2);
xlsx2 = strcat(results.FP.path, '/Data/AUCdata.xlsx');

MaxValue = array2table(MaxValue, 'VariableNames', colnames2);
xlsx3 = strcat(results.FP.path, '/Data/Peaksdata.xlsx');

    writetable(FinalData, xlsx1, 'WriteVariableNames', true); 
    writetable(AUC, xlsx2, 'WriteVariableNames', true); 
    writetable(MaxValue, xlsx3, 'WriteVariableNames', true); 


% Save Figures (Heatmap + DF_Zscore trace):
% Save Figures in same path (figures file should be already created from FFPreprocsingle):

figurepath = strcat(results.FP.path, '/Figures/ZScore_Heatmap.jpg');
saveas(fig1, figurepath);


end % de la función