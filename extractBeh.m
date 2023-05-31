
function [results] = extractBeh(behdatapath, behtype, resultspath, varargin)

% Extract or process behavioral data information and save in same struct

% INPUTS:
%   behdatapath: depending on behavior user must input the folderpath
%   (freezing) or TDTdata filepath obtained after applying fp2mat.m (wheel).
%   behtype: freezing or wheel
%   resultspath: complete filepath of results.mat (struct saved by fp_preproc.m)

% OPTIONAL INPUTS (only used in case of 'wheel' experiment):
%   Movement threshold: set at 2,5 cm/s as a default. 
%   VR: yes or no

% OUTPUT:
%   results: struct with all data

p = inputParser;
addParameter(p,'MovementThreshold', 2.5,@isnumeric);
addParameter(p,'containsVR', "no", @ischar);
parse(p,varargin{:});
MovementThreshold = p.Results.MovementThreshold;
containsVR = p.Results.containsVR;

%% 1. Select behavior 

if behtype == "freezing" | behtype == "Freezing"
    behtype = 1;
elseif behtype == "wheel" | behtype == "Wheel"
    behtype = 2;
end


%% 2.1. Freezing:

if behtype == 1

    BHDpath = behdatapath; % In this case the path is a folder

    % Get path for BHD *.mat files.

    Behpath = strcat(BHDpath, '/Behavior.mat');
    Metpath = strcat(BHDpath, '/Metrics.mat');
    Fspath = strcat(BHDpath, '/Params.mat');

    % Load files*.mat to function workspace

    load(Behpath);
    load(Metpath);
    load(Fspath);


    % Load results file to function workspace. This is where this variables will be stored 
    % (the same file that was created after fp_preproc):

    load(resultspath, 'results')

    % Save variables to struct results.

    results.Behavior.Event.Time = Behavior.Freezing.Bouts;
    results.Behavior.Event.Type = 'Freezing';
    results.Behavior.Location = Metrics.Location.';
    results.Behavior.Fs = Params.Video.frameRate;
    results.Behavior.DistanceTraveled = Metrics.Movement.DistanceTraveled.';
    results.Behavior.Velocity = [Metrics.Movement.Data(:,1), Metrics.Movement.Data].';

    % Save file in same path

    save(resultspath, 'results');

end



%% 2.2. Wheel:

if behtype == 2

    load(behdatapath, 'TDTdata') % In this case behdatapath is a complete filepath

% 2.2.1. Plot and select which is the real wheel signal (it can vary)

    figure(1)
    T = tiledlayout(size(TDTdata.streams.Fi1r.data, 1)+1, 1);

    nexttile
    plot(TDTdata.streams.Wav1.data)
    hold on
    title('Wav1')

    for i = 1:size(TDTdata.streams.Fi1r.data, 1)
        nexttile
        plot(TDTdata.streams.Fi1r.data(i,:))
        hold on
        name = sprintf('Fi1r Row %i', i);
        title(name)
        hold off
    end


% 2.2.2. Select correct wheel channel:

uiwait(msgbox('Select apropiate wheel variable', 'Instructions', "modal")); 

list = {'Wav1','Fi1r Row 1','Fi1r Row 2','Fi1r Row 3','Fi1r Row 4'};
[indx,tf] = listdlg('ListString',list);

if indx == 1 
        wheel = TDTdata.streams.Wav1;
        wvector = TDTdata.streams.Wav1.data;
        wheelfs = TDTdata.streams.Wav1.fs;
elseif indx == 2
        wheel = TDTdata.streams.Fi1r;
        wvector = TDTdata.streams.Fi1r.data(1, :);
        wheelfs = TDTdata.streams.Fi1r.fs;
elseif indx == 3
        wheel = TDTdata.streams.Fi1r;
        wvector = TDTdata.streams.Fi1r.data(2, :);
        wheelfs = TDTdata.streams.Fi1r.fs;
elseif indx == 4
        wheel = TDTdata.streams.Fi1r;
        wvector = TDTdata.streams.Fi1r.data(3, :);
        wheelfs = TDTdata.streams.Fi1r.fs;
elseif indx == 5
        wheel = TDTdata.streams.Fi1r;
        wvector = TDTdata.streams.Fi1r.data(4, :);
        wheelfs = TDTdata.streams.Fi1r.fs;
end

% Load results.mat to function workspace to add behdata:

load(resultspath, 'results')
results.Behavior.Event.Type = 'Wheel';

% 2.2.3. Process wheel signal:

wheeltime = ([0:(length(wheel.data)-1)])./wheelfs; % En segundos 
wheeldist = pi*40;

% Modified JULIO ESPARZA'S CODE:

Lwheel = wheeldist;
speed_movmean_win = 100;

% Filt
posWheel = wvector; 

posWheel = lowpass(posWheel, 1, wheelfs); 
posWheel = (posWheel - min(posWheel)) / (max(posWheel) - min(posWheel)) * Lwheel; 

posWheelSpeed = movmean(posWheel,speed_movmean_win); 
speedIdxIgnore = movmean(diff(posWheel)*wheelfs<-50,round(wheelfs/10))>0; % Detect peaks
speedWheel = abs(diff(movmean(posWheelSpeed,speed_movmean_win)))*wheelfs; % Calculate instant speed
speedWheel(speedIdxIgnore) = nan; 

speedWheel = fillmissing(speedWheel,'spline'); % Interpolate peak values
speedMax = max(diff(movmean(posWheelSpeed(1:end),speed_movmean_win))*wheelfs); 


speedWheel(speedWheel>speedMax) = speedMax;  
speedWheel = [0, speedWheel];  
speedWheel(speedWheel<0) = 0; 
posWheel(posWheel>Lwheel) = Lwheel; 
posWheel(posWheel<0) = 0; 


% 2.2.4. Movement threshold:

threshold = MovementThreshold;

movidx = speedWheel > threshold; 

movidx = find(movidx ~= 0); 
nonadjacent = [0,find(diff(movidx)>1),length(movidx)]; 

% Find events where threshold is met:

for k=2:length(nonadjacent)
    bouts{k-1} = movidx(nonadjacent(k-1)+1):movidx(nonadjacent(k));
    mov_onset(k-1) = min(bouts{k-1}); % Timestamp onset
    mov_offset(k-1) = max(bouts{k-1}); % Timestamp offset
end

movement_bouts = [mov_onset; mov_offset].';


% 2.2.5. Wheel data has to be the same length as FP data (FP data was arranged
% in fp_preproc.m

% Take params data and adapt wheel vector to the same length as FP data: 

SampleBuffer = results.FP.params.SampleBuffer;

preSB = SampleBuffer(1);
postSB = SampleBuffer(2);

pre = round(preSB * wheelfs); 
post = round(postSB * wheelfs);

posWheel = posWheel(pre + 1 : end - post);
speedWheel = speedWheel(pre + 1 : end - post);
wvector = wvector(pre + 1 : end - post);


% 2.2.6. Contains VR?

if containsVR == 'yes'
    hasVR = 1;
elseif containsVR == 'no'
    hasVR = 0;
end


if hasVR == 1
    prompt = {'Start Time of VR (s): '};

    dlgtitle = 'VR Time';
    dims = [1 35];
    definput = {'450'};

    VRtime = inputdlg(prompt, dlgtitle, dims, definput);

end

if hasVR == 1

VRtime = str2num(VRtime{1});

results.Behavior.params.VR.log = 1;
results.Behavior.params.VR.starttimeinsec = VRtime;

    if VRtime > preSB 
        VRtime = VRtime - preSB;
        samplesSignal = round(results.FP.params.fs * VRtime);
        samplesWheel = round(wheelfs * VRtime);
    elseif VRtime <= preSB
        samplesSignal = 0;
        samplesWheel = 0;
    end

% How many samples to cut? 

    if samplesSignal > 0

% Apply to all signals (original like GCAMP and processed like DFF, etc):

            posWheelVR = posWheel(samplesSignal + 1: end); % Posición ya procesada
            speedWheelVR = speedWheel(samplesSignal + 1: end); % Velocidad
            wvectorVR = wvector(samplesSignal + 1: end); % Posición sin procesar con valores originales

% Save VR cut signals:

            results.FP.SignalsVR.DFFModZscore = results.FP.Signals.DFFModZscore(samplesSignal + 1: end);
            results.FP.SignalsVR.DFF = results.FP.Signals.DFF(samplesSignal + 1: end);

            results.FP.SignalsVR.raw.GCAMP = results.FP.Signals.raw.GCAMP(samplesSignal + 1: end);
            results.FP.SignalsVR.raw.Time = results.FP.Signals.raw.Time(samplesSignal + 1: end);
            results.FP.SignalsVR.raw.Isos = results.FP.Signals.raw.Isos(samplesSignal + 1: end);

            results.FP.SignalsVR.smoothed.GCAMP = results.FP.Signals.smoothed.GCAMP(samplesSignal + 1: end);
            results.FP.SignalsVR.smoothed.Time = results.FP.Signals.smoothed.Time(samplesSignal + 1: end);
            results.FP.SignalsVR.smoothed.Isos = results.FP.Signals.smoothed.Isos(samplesSignal + 1: end);

            results.Behavior.PositionVR.original = wvectorVR;
            results.Behavior.PositionVR.processed = posWheelVR;
            results.Behavior.SpeedVR = speedWheelVR;

            results.FP.params.VRstartsample = samplesSignal;
            results.Behavior.params.VRstartsample = samplesWheel;

% Save no VR signals from same recording: 

            posWheelnoVR = posWheel(1:samplesSignal); 
            speedWheelnoVR = speedWheel(1:samplesSignal); 
            wvectornoVR = wvector(1:samplesSignal); 


            results.FP.SignalsnoVR.DFFModZscore = results.FP.Signals.DFFModZscore(1:samplesSignal);
            results.FP.SignalsnoVR.DFF = results.FP.Signals.DFF(1:samplesSignal);

            results.FP.SignalsnoVR.raw.GCAMP = results.FP.Signals.raw.GCAMP(1:samplesSignal);
            results.FP.SignalsnoVR.raw.Time = results.FP.Signals.raw.Time(1:samplesSignal);
            results.FP.SignalsnoVR.raw.Isos = results.FP.Signals.raw.Isos(1:samplesSignal);

            results.FP.SignalsnoVR.smoothed.GCAMP = results.FP.Signals.smoothed.GCAMP(1:samplesSignal);
            results.FP.SignalsnoVR.smoothed.Time = results.FP.Signals.smoothed.Time(1:samplesSignal);
            results.FP.SignalsnoVR.smoothed.Isos = results.FP.Signals.smoothed.Isos(1:samplesSignal);

            results.Behavior.PositionnoVR.original = wvectornoVR;
            results.Behavior.PositionnoVR.processed = posWheelnoVR;
            results.Behavior.SpeednoVR = speedWheelnoVR;



      end

end

% Save variables to struct results. 

results.Behavior.Event.Time = movement_bouts;
results.Behavior.Event.Type = 'Wheel';
results.Behavior.Position.processed = posWheel;
results.Behavior.Speed = speedWheel;
results.Behavior.Position.original = wvector;
results.Behavior.Fs = wheelfs;
results.Behavior.params.MovThreshold = threshold;


% Save file in same path 

save(resultspath, 'results');

end 

end 