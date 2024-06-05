
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% phyPhot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pipeline for processing fiber photometry signal and obtaining descriptive analysis 
% as a result of the combination with behavioral and other physiological 
% signals (Local Field Potential).

% This document contains the script to execute the functions created for
% this purpose and examples of its implementation


%% 1. Load and save TDT data in selected folders 

% You can either write the full path (L17) where TDT data was stored or browse it
% interactively (uigetdir):

% 'mypath' is a character vector of the folder path containting TDT files:

mypath = 'C:\Users\PC\OneDrive - Universidad Complutense de Madrid (UCM)\fotometria\Calb-GCAMP_2AD\191024\Calb-GCAMP_2AD-191024-175908';
filename = 'fpdata.mat';
mypath = uigetdir();
TDTdata = fp2mat(mypath, filename); 

% For freezing files with pulse index (indicating the beginning of fear
% conditioning test), execute the following function:

% TDTdata = cut2wav(); % It cuts the signal specifically to 5 minutes of behavior

% TDTdata is loaded to current workspace

%% 2. Extract, process and plot raw data to detect noise, major artifacts, etc. 

results = fp_preproc(TDTdata, 'SampleBuffer', [5 5]);
results = artifact_fp_preproc(TDTdata, 'SampleBuffer', [5 5]);

%% 3. Extract behavioral data separately:

% Example for wheel behavior:

behdatapath = uigetdir();
behdatapath = strcat(behdatapath, '/', filename);
behtype = 'wheel';
resultspath = uigetdir();
resultspath = strcat(resultspath, '/results.mat');
results = extractBeh(behdatapath, behtype, resultspath); 

% Example for freezing behavior:

behdatapath = uigetdir();
%behdatapath = strcat(behdatapath, '/', filename);
behtype = 'freezing';
resultspath = uigetdir();
resultspath = strcat(resultspath, '/results.mat');
results = extractBeh(behdatapath, behtype, resultspath);


% for ii = 1:size(TrialData, 1)
%     plot((1:51), TrialData(ii, :), 'blue')
%     hold on
% end
% 
% eventsfs = results.Behavior.Event.Time./BehFs;
% eventsfs = round(eventsfs.*Fs);
% 
% for ii = 1:size(eventsfs, 1)-1
%     plot(DFFZ(eventsfs(ii, 1):(eventsfs(ii, 2))), 'blue')
%     hold on
% end
% 
% plot(CSBL)
%% 4. Combine behavioral and FP data:

% 4.1. Freezing:

% Plot PETH and calculate AUC and Peaks:

results = PETHZscore(resultspath, 'baseline', 5, 5, 0.2);
% Fs =results.FP.params.fs;
% Ts = results.FP.Signals.raw.Time;
% Ch490 = results.FP.Signals.raw.GCAMP;
% Ch405 = results.FP.Signals.raw.Isos;
% Pre = 2;
% Post = 2;
% bin = 0.2;
% BehFs = 14.99;
% EventTS = results.Behavior.Event.Time./BehFs;
% BaselineWind = 5;
% BaselineWind2 = 5;
% which TrialPETHZScore
% FigTitle = '';
% Bin = 30;
% BL_Width = 2;
% BL_Start = 2;
% TrialPETHZScore(Fs, Ts, Pre, Post, EventTS, Ch490,Ch405,'AUC','',5,4,5);
%         TrialPETHZScore(Fs, Ts, Pre,...
%             Post, EventTS, Ch490, Ch405,...
%             '',FigTitle,Bin,0,5, BL_Width, BL_Start);
% 
% 4.2. Wheel:

% Descriptives of wheel position and FP:

results = fpbinpos(resultspath);

% Descriptives of wheel position, FP and LFP:

fppath = strcat(mypath, '/', filename);
results = FP_LFP(fppath, resultspath);



% 5. ONSET Y OFFSET DE LOS EVENTOS

% con valores predefinidos
results = dffPETHonset(results);
results = dffPETHoffset(results);

% ajustando valores :

results = dffPETHonset(results, ...
    'Pre', 3, ... % ventana de tiempo pre en segundos
    'Post', 3, ... % ventana de tiempo post en segundos
    'bin', 0.2, ... % la fs del PETH
    'binthreshold', 1, ... % cuántos seg tiene que durar como mínimo el baseline
    'AUCqt', 3, ... % cantidad de intervalos a evaluar para sacar AUC y density
    'AUCint', [-2 -1; 0 1; 1 2], ... % períodos de intervalos de principio a fin. Si es pre, en negativo. % Tiene que ser una matriz donde cada fila es un intervalo y cada columna es el principio y el fin del intervalo (2 cols)
    'savepath', yourpath, ... % Por si lo quieres guardar en otro sitio
    'summarymeasures', ["median" "sem"]); % Por si quieres utilizar otras medidas resumen de los datos

% lo mismo para offset, teniendo en cuenta que en este caso:
% >> Pre: es el episodio de freezing
% >> Post: es el baseline

% Default params: 

% addParameter(p,'Pre', 2,@isnumeric);
% addParameter(p,'Post', 2, @isnumeric);
% addParameter(p,'bin', 0.1, @isnumeric);
% % addParameter(p,'BinWinStart', 2, @isnumeric);
% % addParameter(p,'BinWinDuration', 1.5, @isnumeric);
% addParameter(p,'binthreshold', 0.5, @isnumeric);
% addParameter(p,'AUCqt', 2, @isnumeric);
% addParameter(p,'AUCint', [-2 0; 0 2], @ismatrix);
% addParameter(p,'savepath', results.FP.path, @ischaracter);
% addParameter(p,'summarymeasures', ["mean", "std"], @isstring);





