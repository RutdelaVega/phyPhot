
function [results] = dffPETHoffset(results, varargin)

%% dffPETHoffset plots and analyzes end of behavior events

p = inputParser;
addParameter(p,'Pre', 2,@isnumeric);
addParameter(p,'Post', 2, @isnumeric);
addParameter(p,'bin', 0.1, @isnumeric);
addParameter(p,'BLthreshold', 0.5, @isnumeric);
addParameter(p,'AUCqt', 2, @isnumeric);
addParameter(p,'AUCint', [-2 0; 0 2], @ismatrix);
addParameter(p,'savepath', results.FP.path, @ischaracter);
addParameter(p,'summarymeasures', ["mean", "std"], @isstring);

parse(p,varargin{:});
Pre = p.Results.Pre;
Post = p.Results.Post;
bin = p.Results.bin;
threshold = p.Results.BLthreshold;
AUCqt = p.Results.AUCqt;
AUCint = p.Results.AUCint;
savepath = p.Results.savepath;
summarymeasures = p.Results.summarymeasures;

%% 1. Load data: 


Fs = results.FP.params.fs;
BehFs = results.Behavior.Fs;
evsamp = results.Behavior.Event.Time;
Events = evsamp./BehFs; 
if Events(end, 2) > 300
    Events(end, 2) =  300; % en seg
    evsamp(end, 2) = 4500; % en BehFs
end

DFFZ = double(results.FP.Signals.DFFModZscore);
DFF = results.FP.Signals.DFF;
Time = results.FP.Signals.raw.Time;

%% 2. Variables for PETH: 

 binsize = round(bin*Fs); % Unit is in samples, if each bin represents 0.2 seconds, then how many samples are within 0.2 seconds (in FP fs)
 PreWind = round(Pre*Fs);
 PostWind = round(Post*Fs);
% Transform timestamps in BHD to timestamps in GCAMP:
FPtimestamps = round(Events.*Fs);



%% 3. Save which events have overlapping windows

evdif = Events(2:end, 1) - Events(1:end-1, 2); % Diferencia entre eventos (idx +1)
overlap =  evdif < Post;
boutdurTS = FPtimestamps(:, 2) - FPtimestamps(:, 1);



%% 4. Max length of freezing event >> esto es la PreWIND (ahora está antes que el BL): 

boutdur = Events(:, 2) - Events(:, 1);
[val, ~] = max(boutdur); Prewfr = val; PreWind = round(Prewfr*Fs);

EVdata = zeros([size(Events, 1) PreWind]);
BLdata = zeros([size(Events, 1) PostWind]);


%% 5. Extract FP data from each event and save data:

EVDFF = EVdata;
BLDFF = BLdata;

nansizepre = abs(boutdurTS - size(EVdata, 2)) - 1;


% COMPLETE DATA WITH OVERLAP

for ii = 1:length(Events)
    if ii == 1
        BLDFF(ii, :) = DFF(FPtimestamps(ii, 2)+1:FPtimestamps(ii, 2)+PostWind);
        if (FPtimestamps(ii, 1) - PreWind) < 0 % que exceda por deltante
            % añadir tantos nans por delante == desde inicio de freezing hasta
            % distancia que haya para llegar a PreWind muestras (osea PreWind -
            % las muestras que haya en ese intervalo >> boutdurTS) PreWind -
            % boutdurTS
            EVDFF(ii, :) = [nan([1 nansizepre(ii)]) DFF(FPtimestamps(ii, 1):FPtimestamps(ii, 2))];
        else % caso normal
            EVDFF(ii, :) = DFF(FPtimestamps(ii, 2)-PreWind+1:FPtimestamps(ii, 2));
        end
    elseif ii >1 && ii < length(Events)
        if (FPtimestamps(ii, 1) - PreWind) < 0 % Que exceda por delante
            BLDFF(ii, :) = DFF(FPtimestamps(ii, 2)+1:FPtimestamps(ii, 2)+PostWind);
            EVDFF(ii, :) = [nan([1 nansizepre(ii)]) DFF(FPtimestamps(ii, 1):FPtimestamps(ii, 2))];
        elseif (FPtimestamps(ii, 2) +1 + PostWind) > length(DFF) % Que exceda por detrás
            nansize = (PostWind - length(DFF(FPtimestamps(ii, 2)+1:end)));
            BLDFF(ii, :) = [DFF(FPtimestamps(ii, 2)+1:end) nan([1 nansize])];
            EVDFF(ii, :) = DFF(FPtimestamps(ii, 2)-PreWind+1:FPtimestamps(ii, 2));
            % añadir nans por el final
        else % caso normal
            BLDFF(ii, :) = DFF(FPtimestamps(ii, 2)+1:FPtimestamps(ii, 2)+PostWind);
            EVDFF(ii, :) = DFF(FPtimestamps(ii, 2)-PreWind+1:FPtimestamps(ii, 2));
        end

    elseif ii == length(Events)
        if (FPtimestamps(ii, 2) +1 + PostWind) > length(DFF)
            nansize = (PostWind - length(DFF(FPtimestamps(ii, 2)+1:end)));
            BLDFF(ii, :) = [DFF(FPtimestamps(ii, 2)+1:end) nan([1 nansize])];
            EVDFF(ii, :) = DFF(FPtimestamps(ii, 2)-PreWind+1:FPtimestamps(ii, 2));            
        else % caso normal
            BLDFF(ii, :) = DFF(FPtimestamps(ii, 2)+1:FPtimestamps(ii, 2)+PostWind);
            EVDFF(ii, :) = DFF(FPtimestamps(ii, 2)-PreWind+1:FPtimestamps(ii, 2));
        end
    end
end


%% 6. Correct for length: assign nans where there is no freezing

% BASELINE: put nans for overlapping trials: (first event is already
% corrected)

overlapsize = zeros([size(overlap) 1]);
overlapsize(overlap) = evdif(overlap) - Post;
overlapsize = [overlapsize; 0];
overlapsize = abs(round(overlapsize.*Fs));

nanBLDFF = BLDFF;


for jj = 1:size(BLDFF, 1)
    if overlapsize(jj) == 0
        continue
    else
    nanvect = nan([1 overlapsize(jj)]);
    nanBLDFF(jj, PostWind-overlapsize(jj)+1:PostWind) = nanvect;
    end
end


% EVENT:

nanEVDFF = EVDFF; % Estas vars ahora contienen nan, de ahí el nombre

boutdursamp = round(boutdur.*Fs);
elimsample = size(EVDFF, 2) - boutdursamp;

for kk = 1:size(EVDFF, 1) 
    nanvect = nan([1 elimsample(kk)]);
    nanEVDFF(kk, 1:size(EVDFF, 2)-boutdursamp(kk)) = nanvect;
end




%% 7. Compute DFF

for nn = 1:size(BLdata)

    BLmed(nn) = median(nanBLDFF(nn, :), 'omitnan'); BLmad(nn) = mad(nanBLDFF(nn, :));
    EvDFFZ(nn, :) = (nanEVDFF(nn, :) - BLmed(nn))./BLmad(nn);  % DFFZ on baseline
% WITH NORMALIZATION OF BL:
     BLDFFZ(nn, :) = (nanBLDFF(nn, :) - BLmed(nn))./BLmad(nn);
end


%% 8. Downsample through binsize for PETH:

BLjumps = 1:binsize:size(BLDFFZ, 2);
EVjumps = 1:binsize:size(EvDFFZ, 2);

for ii = 1:length(BLjumps)-1
    binnedBL(:,ii) = median(BLDFFZ(:, BLjumps(ii):BLjumps(ii+1)), 2, 'omitnan');
end


for ii = 1:length(EVjumps)-1
    binnedEV(:, ii) = median(EvDFFZ(:, EVjumps(ii):EVjumps(ii+1)), 2, 'omitnan');
end

%% 9. Prepare Data for PETH:


% Criteria for BL (threshold = 0.5 seg):

lenBL = zeros([size(BLDFFZ, 1), 1]);
for kk = 1:size(BLDFFZ, 1)
    lenBL(kk) = sum(~isnan(BLDFFZ(kk, :)));
end

sampthresh = round(threshold*size(BLDFFZ, 2)./Post);
underthresh = lenBL < sampthresh;

% Data and plot PETH:

PETHdata = [binnedEV binnedBL];

PETHfreez = PETHdata; % Save to other var to be saved later wo the criteria and full data

PETHdata = PETHdata(~underthresh, :);
PETHts = linspace(-Prewfr, Post, size(PETHdata, 2)); 
PETHevnr = size(PETHdata, 1);


% Extract summary measures:

if summarymeasures(1) == "mean"
    PETHm = mean(PETHfreez, 1, 'omitnan');
elseif summarymeasures(1) == "median"
    PETHm = median(PETHfreez, 1, 'omitnan');
end

if summarymeasures(2) == "std"
    PETHd = std(PETHfreez, 1, 'omitnan');
elseif summarymeasures(2) == "sem"
    PETHd = std(PETHfreez, 1, 'omitnan')./sqrt(size(PETHfreez, 1));
elseif summarymeasures(2) == "ci"
    SEM = td(PETHfreez, 1, 'omitnan')./sqrt(size(PETHfreez, 1));               % Standard Error
    ts = tinv([0.025  0.975],size(PETHfreez, 1)-1);      % T-Score
    PETHd = mean(PETHfreez, 1, 'omitnan') + ts*SEM;                      % Confidence Intervals
elseif summarymeasures(2) == "mad"
    PETHd = mad(PETHfreez);
end




fig1 = figure(1);
subplot(2, 1, 1)
h = imagesc(PETHts, 1:PETHevnr, (PETHdata));
set(h, 'AlphaData', ~isnan(PETHdata))
clim([-5 5])
subplot(2, 1, 2)
PETHmedia = mean(PETHdata, 1, 'omitnan');
g = plot(PETHts, median(PETHdata, 1, 'omitnan'), 'LineWidth', 2);
xlim([min(PETHts) max(PETHts)])
xline(0, '-r')
yline(mean(mean(BLDFFZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)


fig2 = figure(2);
plot(Time, DFFZ)
hold on
for ii = 1:size(Events, 1)
    xline(Events(ii, 1), '-g')
    xline(Events(ii, 2), '-r')
end
title('ZScore DFF with onset freezing events')


% Sort freezing episodes by length: 

[~, idx] = sort(boutdur, 'ascend');
%in case of repeated durations:

for  ii= 1:size(boutdur, 1)
    sortedPETH(ii, :) = PETHfreez((idx(ii)), :);
end


% Criteria for BL (threshold = 0.5 seg):

rowstoelim = underthresh == 1;
if ~isempty(rowstoelim)
    sortedPETH = sortedPETH(~rowstoelim, :);
end

fig3 = figure(3);
tiledlayout(2, 1)
nexttile
h = imagesc(PETHts, 1:PETHevnr, (sortedPETH));
set(h, 'AlphaData', ~isnan(sortedPETH))
clim([-5 5])
nexttile
plot(PETHts, median(PETHdata, 1, 'omitnan'), 'LineWidth', 2);
xlim([min(PETHts) max(PETHts)])
xline(0, '-r')
yline(mean(mean(BLDFFZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)

% Recortar PETH a rangos deseados: 

idxPETHts = PETHts >= -Pre & PETHts <= Post;
udPETHts = PETHts(idxPETHts);
udPETHdata = PETHdata(:, idxPETHts);

fig4 = figure(4);
tiledlayout(2, 1)
nexttile
h = imagesc(udPETHts, 1:size(udPETHdata, 1), udPETHdata);
set(h, 'AlphaData', ~isnan(udPETHdata))
clim([-5 5])
nexttile
plot(udPETHts, median(udPETHdata, 1, 'omitnan'), 'LineWidth', 2);
xlim([min(udPETHts) max(udPETHts)])
%ylim([min(udPETHdata) max(udPETHdata)])
xline(0, '-r')
yline(mean(mean(BLDFFZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)


if summarymeasures(1) == "mean"
    udPETHm = mean(udPETHdata, 1, 'omitnan');
elseif summarymeasures(1) == "median"
    udPETHm = median(udPETHdata, 1, 'omitnan');
end

if summarymeasures(2) == "std"
    udPETHd = std(udPETHdata, 1, 'omitnan');
elseif summarymeasures(2) == "sem"
    udPETHd = std(udPETHdata, 1, 'omitnan')./sqrt(size(udPETHdata, 1));
elseif summarymeasures(2) == "ci"
    SEM = td(udPETHdata, 1, 'omitnan')./sqrt(size(udPETHdata, 1));               % Standard Error
    ts = tinv([0.025  0.975],size(udPETHdata, 1)-1);      % T-Score
    udPETHd = mean(udPETHdata, 1, 'omitnan') + ts*SEM;                      % Confidence Intervals
elseif summarymeasures(2) == "mad"
    udPETHd = mad(udPETHdata);
end


%% 11. AUC + Peak + density of AUC: 

PETHfs = bin;
PETHtsAUC = PETHts + Prewfr; % Escalo a 0

for jj = 1:size(binnedBL, 1)
    if double(sum(~isnan(binnedBL(jj, :))))*PETHfs < 0.5 % criterio de duración.. si ponemos una fs para el PETH muy alto, hay que tener cuidado porque podría pasar mucho esto al ser un ds
        tbl(jj) = sum(~isnan(binnedBL(jj, :)))*PETHfs; AUCbl(jj) = nan;
        densbl(jj) = nan; AUCev(jj) = nan; densev(jj) = nan;
    else
    tbl(jj) = double(sum(~isnan(binnedBL(jj, :))))*PETHfs; % duración bin en seg // la dur del bout ya la tenemos
    tempBLts = PETHtsAUC(1:sum(~isnan(binnedBL(jj,:))));
    AUCbl(jj) = trapz(tempBLts, binnedBL(jj, ~isnan(binnedBL(jj, :))));
    densbl(jj) = AUCbl(jj)./length(binnedBL(jj, ~isnan(binnedBL(jj, :))));
    tempEVts = PETHtsAUC(1:sum(~isnan(binnedEV(jj, :))));
    AUCev(jj) = trapz(tempEVts,binnedEV(jj, ~isnan(binnedEV(jj, :))));
    densev(jj) = AUCev(jj)./length(binnedEV(jj, ~isnan(binnedEV(jj, :)))); 
    end    
end


% 11.2. Por intervalos de tiempo determinados por usuario: 

AUCdata = zeros([length(boutdur) size(AUCint, 1)]);
DENSdata = AUCdata;
AUCdur = AUCdata;

for ii = 1:AUCqt
    for jj = 1:size(binnedBL, 1)
        if double(sum(PETHts >= AUCint(ii, 1) & PETHts < AUCint(ii, 2) & ~isnan(PETHfreez(jj, :))))*bin < 0.5
            AUCdata(jj, ii) = nan;
            DENSdata(jj, ii) = nan;
            AUCdur(jj, ii) = double(sum(PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHfreez(jj, :))))*bin;
        else
        AUCdata(jj, ii) = trapz(bin, PETHfreez(jj, (PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHfreez(jj, :)))));
        AUCdur(jj, ii) = double(sum(PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHfreez(jj, :))))*bin;
        DENSdata(jj, ii) = AUCdata(jj, ii)./length(PETHfreez(jj, (PETHts >= AUCint(ii, 1) & PETHts < AUCint(ii, 2) & ~isnan(PETHfreez(jj, :)))));
        end
    end
end  

% 11. 3. Filtrar los que cumplan la condición que todos tengan la duración
% indicada por usuario. 

intdur = repmat((AUCint(:, 2) - AUCint(:, 1)).', length(boutdur), 1);
indx = double(AUCdur == intdur); 
indx = sum(indx, 2);
AUCsamelen = AUCdata(indx == AUCqt, :); % keep only those where all events are the chosen length by user (not necessarily the same length between intervals)
% It could happen that user chooses 1int of 2 s and another of 3 s. 
DENSsamelen = DENSdata(indx == AUCqt, :);
AUCdursamelen = AUCdur(indx == AUCqt, :); 
boutdursamelen = boutdur(indx == AUCqt, :);

% Measurenames: 

% Data 1:

measurenames1 = {'AUCpre', 'AUCpost', 'DENSpre (Apre/n)', 'DENSpost(Apost/n)', 'durpre (s)', 'boutdur (s)'};

prepostFP = table(AUCbl.', AUCev.', densbl.', densev.', tbl.', boutdur, 'VariableNames', measurenames1);


% Data 2 and 3 (even if empty):

AUCnames = {};
DENSnames = {};
durnames = {};

for ii = 1:AUCqt % 3 because of three times
    AUCnames{ii} = sprintf('AUC Int%i', ii);
    DENSnames{ii} = sprintf('DENS Int%i', ii);
    durnames{ii} = sprintf('Dur Int%i', ii);
    
end

% if ~isempty(AUCsamelen)
%     for ii = 1:size(AUCsamelen)
%         AUCsamenames{ii} = sprintf('AUC Int%i', ii);
%         DENSsamenames{ii} = sprintf('DENS Int%i', ii)
%     end
% end
% 

measurenames2 = [AUCnames DENSnames durnames 'boutdur'];
allintFP = array2table([AUCdata DENSdata AUCdur boutdur], 'VariableNames', measurenames2);
if ~isempty(AUCsamelen)
    samelenFP = array2table([AUCsamelen DENSsamelen AUCdursamelen boutdursamelen], 'VariableNames', measurenames2);
    %xlsx3 = writetable(samelen);
end

%% 12. Save all data and figures: 

% Params:
results.PETH.offset.params.Pre = Pre;
results.PETH.offset.params.Post = Post;
results.PETH.offset.params.Prewfr = Prewfr;
results.PETH.offset.params.PreWind = PreWind;
results.PETH.offset.params.PostWind = PostWind;
results.PETH.offset.params.bin = bin;
results.PETH.offset.params.binsize = binsize; % in samples
results.PETH.offset.params.binthreshold = threshold;
results.PETH.offset.params.summarymeasures = summarymeasures;

% PETH variables: 

results.PETH.offset.PETHdata = PETHfreez;
results.PETH.offset.PETHwbincriteria = PETHdata;
results.PETH.offset.PETHts = PETHts;
results.PETH.offset.udPETHdata = udPETHdata;
results.PETH.offset.udPETHts = udPETHts;

results.PETH.offset.PETHmean = PETHm;
results.PETH.offset.PETHdev = PETHd;
results.PETH.offset.udPETHmean = udPETHm;
results.PETH.offset.udPETHdev = udPETHd;

results.PETH.offset.sortedPETH = sortedPETH;
results.PETH.offset.EvDFFZ = EvDFFZ;
results.PETH.offset.BLDFFZ = BLDFFZ;



% Figures: 

folderpath1 = strcat(savepath, '/Figures');
folderpath2 = strcat(savepath, '/Data');
mkdir(folderpath1, 'offsetPETH');
mkdir(folderpath2, 'offsetPETH')

fig1name = strcat(results.FP.path, '/Figures/offsetPETH/offcompletePETH.jpeg');
fig2name = strcat(results.FP.path, '/Figures/offsetPETH/offeventsDFF.jpeg');
fig3name = strcat(results.FP.path, '/Figures/offsetPETH/offsortedPETH.jpeg');
fig4name = strcat(results.FP.path, '/Figures/offsetPETH/offblthreshPETH.jpeg');

saveas(fig1, fig1name)
saveas(fig2, fig2name)
saveas(fig3, fig3name)
saveas(fig4, fig4name)

% Tables:

xlsx1 = strcat(results.FP.path, '/Data/offsetPETH/offprepostAUC.xlsx');
xlsx2 = strcat(results.FP.path, '/Data/offsetPETH/offallintAUC.xlsx');
xlsx3 = strcat(results.FP.path, '/Data/offsetPETH/offsamelengthintAUC.xlsx');
writetable(prepostFP, xlsx1, 'WriteVariableNames', true); 
writetable(allintFP, xlsx2, 'WriteVariableNames', true);
if ~isempty(AUCsamelen)
    writetable(samelenFP, xlsx3, 'WriteVariableNames', true);
% else
%     message = 'Not enough data from chosen interval (data is not the same length). Data Matrix is empty.';
%     disp(message)
end

% Events:

pethcolcompletets = string(PETHts);
xlsx4 = strcat(results.FP.path, '/Data/offsetPETH/completePETH.xlsx');
xlsx5 = strcat(results.FP.path, '/Data/offsetPETH/sortedPETHwthresh.xlsx');
xlsx6 = strcat(results.FP.path, '/Data/offsetPETH/PETHwthresh.xlsx');
PETHfreez2 = array2table(PETHfreez, 'VariableNames', pethcolcompletets);
sortedPETH2 = array2table(sortedPETH, 'VariableNames', pethcolcompletets);
PETHdata2 = array2table(PETHdata, 'VariableNames',pethcolcompletets);

writetable(PETHfreez2, xlsx4, 'WriteVariableNames', true);
writetable(sortedPETH2, xlsx5, 'WriteVariableNames', true);
writetable(PETHdata2, xlsx6, 'WriteVariableNames', true);

pethcolintts = string(udPETHts);
xlsx7 = strcat(results.FP.path, '/Data/offsetPETH/udPETH.xlsx');
udPETHdata2 = array2table(udPETHdata, 'VariableNames', pethcolintts);
writetable(udPETHdata2, xlsx7, 'WriteVariableNames', true);



% Mean + sem:

summaryPETH = table(PETHts.', PETHm.', PETHd.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2))});
udsummaryPETH = table(udPETHts, udPETHm, udPETHd, 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2))});

xlsx8 = strcat(results.FP.path, '/Data/offsetPETH/offsummaryPETH.xlsx');
xlsx9 = strcat(results.FP.path, '/Data/offsetPETH/offudsummaryPETH.xlsx');

writetable(summaryPETH, xlsx8)
writetable(udsummaryPETH, xlsx9)




end