
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
mypath = uigetdir()
TDTdata = fp2mat(mypath, filename); 

% For freezing files with pulse index (indicating the beginning of fear
% conditioning test), execute the following function:

TDTdata = cut2wav(); % It cuts the signal specifically to 5 minutes of behavior

% TDTdata is loaded to current workspace

%% 2. Extract, process and plot raw data to detect noise, major artifacts, etc. 

results = fp_preproc(TDTdata);

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

%% 4. Combine behavioral and FP data:

% 4.1. Freezing:

% Plot PETH and calculate AUC and Peaks:

results = PETHZscore(resultspath, 'DFF', 2, 2, 0.2);

% 4.2. Wheel:

% Descriptives of wheel position and FP:

results = fpbinpos(resultspath);

% Descriptives of wheel position, FP and LFP:

fppath = strcat(mypath, '/', filename);
results = FP_LFP(fppath, resultspath);







