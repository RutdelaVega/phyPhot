
function [results] = dffPETHoffset(contvar, fscontvar, evsamp, fsev, combiname, savepath, varargin)

%% dffPETHoffset plots and analyzes end of behavior events starting from the end of the event
% INPUTS:
%   contvar (vector): continuous variable that you want to allign with specific
%  events (z values). DFF, velocity, acceleration, etc. 
%   fscontvar (escalar): frequency sampling of contvar
%   evsamp (two col matrix): matrix of two columns containing onset sample and offset sample of such
%  events. This could be: freezing, shocks, movement, etc. 
%   fsev (escalar): frequency sampling of events, needed as the fs could be different
%  for those two variables. 
%   savepath (string): where to save all the data

% OPTIONAL: 
%  Pre (s): start of PETH window. Ex: 2 is -2 seconds before event. 
%  Post (s): end of PETH window. 2 is duration of event 2 seconds after
%  start
%  bin: size of each square in the PETH
%  threshold (s): criteria to eliminate events (or not show them)
%  AUCqt (n): how many AUC to calculate
%  AUCint (s): intervals to get information from (coherent with PETH
%  window size)
%  summarymeasures (string): which measures to extract for later
%  quantitative analysis.
%  maxlength (s): max length of your recording/file/etc.  


p = inputParser;
addParameter(p,'Pre', 2,@isnumeric);
addParameter(p,'Post', 2, @isnumeric);
addParameter(p,'bin', 0.1, @isnumeric);
addParameter(p,'BLthreshold', 0.5, @isnumeric);
addParameter(p,'AUCqt', 2, @isnumeric);
addParameter(p,'AUCint', [-2 0; 0 2], @ismatrix);
% addParameter(p,'savepath', results.FP.path, @ischaracter);
addParameter(p,'summarymeasures', ["mean", "std"], @isstring);
addParameter(p,'maxlength', 300 ,@isnumeric);

parse(p,varargin{:});
Pre = p.Results.Pre;
Post = p.Results.Post;
bin = p.Results.bin;
threshold = p.Results.BLthreshold;
AUCqt = p.Results.AUCqt;
AUCint = p.Results.AUCint;
summarymeasures = p.Results.summarymeasures;
maxlength = p.Results.maxlength;

%% 1. Load data: 

% 
% Fs = results.FP.params.fs;
% BehFs = results.Behavior.Fs;
% evsamp = results.Behavior.Event.Time;
Events = evsamp./fsev; 
tscontvar = round(Events.*fscontvar);
Time = (0:length(contvar)-1)./fscontvar;
if Events(end, 2) > maxlength
    Events(end, 2) =  maxlength; % en seg
    tscontvar(end, 2) = length(contvar); % en fs de contvar (DFF por ej)
end

if size(contvar, 1) > size(contvar, 2)
    contvar = contvar.';
end

% DFFZ = double(results.FP.Signals.DFFModZscore);
% DFF = results.FP.Signals.DFF;
% Time = results.FP.Signals.raw.Time;

%% 2. Variables for PETH: 

%  binsize = round(bin*Fs); % Unit is in samples, if each bin represents 0.2 seconds, then how many samples are within 0.2 seconds (in FP fs)
%  PreWind = round(Pre*Fs);
%  PostWind = round(Post*Fs);

 binsize = bin*fsev; 
 PreWind = Pre*fsev;
 PostWind = Post*fsev;

 Precontvar = round(Pre*fscontvar);
 Postcontvar = round(Post*fscontvar);

% Transform timestamps in BHD to timestamps in GCAMP:
% FPtimestamps = round(Events.*Fs); % ya está en tscontvar



%% 3. Save which events have overlapping windows

evdif = Events(2:end, 1) - Events(1:end-1, 2); % Diferencia entre eventos (idx +1)
overlap =  evdif < Post;
boutdurTS = tscontvar(:, 2) - tscontvar(:, 1);



%% 4. Max length of freezing event >> esto es la PreWIND (ahora está antes que el BL): 

boutdur = Events(:, 2) - Events(:, 1);
[val, ~] = max(boutdur); Prewfr = val; 
% PreWind = round(Prewfr*Fs);

EVdata = zeros([size(Events, 1) Precontvar]);
BLdata = zeros([size(Events, 1) Postcontvar]);


%% 5. Extract FP data from each event and save data:

nansizepre = abs(boutdurTS - size(EVdata, 2)) - 1;


% COMPLETE DATA WITH OVERLAP

for ii = 1:size(Events, 1)
    if ii == 1
        BLdata(ii, :) = contvar(tscontvar(ii, 2)+1:tscontvar(ii, 2)+Postcontvar);
        if (contvar(ii, 2) - Precontvar) < 0 % que exceda por deltante
            % añadir tantos nans por delante == desde inicio de freezing hasta
            % distancia que haya para llegar a PreWind muestras (osea PreWind -
            % las muestras que haya en ese intervalo >> boutdurTS) PreWind -
            % boutdurTS
            EVdata(ii, :) = [nan([1 nansizepre(ii)]) contvar(tscontvar(ii, 1):tscontvar(ii, 2))];
        else % caso normal
            EVdata(ii, :) = contvar(tscontvar(ii, 2)-Precontvar+1:tscontvar(ii, 2));
        end
    elseif ii >1 && ii < size(Events, 1)
        if (tscontvar(ii, 2) - Precontvar) < 0 % Que exceda por delante
            BLdata(ii, :) = contvar(tscontvar(ii, 2)+1:tscontvar(ii, 2)+Postcontvar);
            EVdata(ii, :) = [nan([1 nansizepre(ii)]) DFF(tscontvar(ii, 1):tscontvar(ii, 2))];
        elseif (tscontvar(ii, 2) +1 + Postcontvar) > length(contvar) % Que exceda por detrás
            nansize = (Postcontvar - length(contvar(tscontvar(ii, 2)+1:end)));
            BLdata(ii, :) = [contvar(tscontvar(ii, 2)+1:end) nan([1 nansize])];
            EVdata(ii, :) = contvar(tscontvar(ii, 2)-Precontvar+1:tscontvar(ii, 2));
            % añadir nans por el final
        else % caso normal
            BLdata(ii, :) = contvar(tscontvar(ii, 2)+1:tscontvar(ii, 2)+Postcontvar);
            EVdata(ii, :) = contvar(tscontvar(ii, 2)-Precontvar+1:tscontvar(ii, 2));
        end

    elseif ii == size(Events, 1)
        if (tscontvar(ii, 2) + Postcontvar) > length(contvar)
            nansize = (Postcontvar - length(DFF(tscontvar(ii, 2)+1:end)));
            BLdata(ii, :) = [contvar(tscontvar(ii, 2)+1:end) nan([1 nansize])];
            EVdata(ii, :) = contvar(tscontvar(ii, 2)-Precontvar+1:tscontvar(ii, 2));            
        else % caso normal
            BLdata(ii, :) = contvar(tscontvar(ii, 2)+1:tscontvar(ii, 2)+Postcontvar);
            EVdata(ii, :) = contvar(tscontvar(ii, 2)-Precontvar+1:tscontvar(ii, 2));
        end
    end
end


%% 6. Correct for length: assign nans where there is no freezing

% BASELINE: put nans for overlapping trials: (first event is already
% corrected)

overlapsize = zeros([size(overlap) 1]);
overlapsize(overlap) = evdif(overlap) - Post;
overlapsize = [overlapsize; 0];
overlapsize = abs(round(overlapsize.*fscontvar));

nanBLdata = BLdata;


for jj = 1:size(BLdata, 1)
    if overlapsize(jj) == 0
        continue
    else
    nanvect = nan([1 overlapsize(jj)]);
    nanBLdata(jj, Postcontvar-overlapsize(jj)+1:Postcontvar) = nanvect;
    end
end


% EVENT:

nanEVdata = EVdata; % Estas vars ahora contienen nan, de ahí el nombre

boutdursamp = round(boutdur.*fscontvar);
elimsample = size(EVdata, 2) - boutdursamp;

for kk = 1:size(EVdata, 1) 
    nanvect = nan([1 elimsample(kk)]);
    nanEVdata(kk, 1:size(EVdata, 2)-boutdursamp(kk)) = nanvect;
end




%% 7. Compute continuous variable

for nn = 1:size(BLdata)

    if sum(~isnan(nanBLdata(nn, :))) < 2
        EVdata(nn, :) = nan([1 size(nanEVdata, 2)]);
        BLdata(nn, :) = nan([1 size(nanBLdata, 2)]);
    else
    BLmed(nn) = median(nanBLdata(nn, :), 'omitnan'); BLmad(nn) = mad(nanBLdata(nn, :));
    EVdataZ(nn, :) = (nanEVDFF(nn, :) - BLmed(nn))./BLmad(nn);  % DFFZ on baseline
% WITH NORMALIZATION OF BL:
    BLdataZ(nn, :) = (nanBLdata(nn, :) - BLmed(nn))./BLmad(nn);
    end
end


%% 8. Downsample through binsize for PETH:

BLjumps = 1:bin*fscontvar:size(BLdataZ, 2);
BLjumps(end+1) = size(BLdataZ, 2);
BLjumps = round(BLjumps);

EVjumps = 1:bin*fscontvar:size(EVdataZ, 2);
EVjumps(end+1) = size(EVdataZ, 2);
EVjumps = round(EVjumps);

for ii = 1:length(BLjumps)-1
    binnedBL(:,ii) = median(BLdataZ(:, BLjumps(ii):BLjumps(ii+1)), 2, 'omitnan');
end

for ii = 1:length(EVjumps)-1
    binnedEV(:, ii) = median(EVdataZ(:, EVjumps(ii):EVjumps(ii+1)), 2, 'omitnan');
end

%% 9. Prepare Data for PETH:


% Criteria for BL (threshold = 0.5 seg):

lenBL = zeros([size(BLdataZ, 1), 1]);
for kk = 1:size(BLdataZ, 1)
    lenBL(kk) = sum(~isnan(BLdataZ(kk, :)));
end

sampthresh = round(threshold*size(BLdataZ, 2)./Post);
underthresh = lenBL < sampthresh;

% Data and plot PETH:

PETHdata = [binnedEV binnedBL];

PETHev = PETHdata; % Save to other var to be saved later wo the criteria and full data

PETHdata = PETHdata(~underthresh, :);
PETHts = linspace(-Prewfr, Post, size(PETHdata, 2)); 
PETHevnr = size(PETHdata, 1);


% Extract summary measures:

if summarymeasures(1) == "mean"
    PETHm = mean(PETHev, 1, 'omitnan');
elseif summarymeasures(1) == "median"
    PETHm = median(PETHev, 1, 'omitnan');
end

if summarymeasures(2) == "std"
    PETHd = std(PETHev, 1, 'omitnan');
elseif summarymeasures(2) == "sem"
    PETHd = std(PETHev, 1, 'omitnan')./sqrt(size(PETHev, 1));
elseif summarymeasures(2) == "ci"
    SEM = td(PETHev, 1, 'omitnan')./sqrt(size(PETHev, 1));               % Standard Error
    ts = tinv([0.025  0.975],size(PETHev, 1)-1);      % T-Score
    PETHd = mean(PETHev, 1, 'omitnan') + ts*SEM;                      % Confidence Intervals
elseif summarymeasures(2) == "mad"
    PETHd = mad(PETHev);
end




fig1 = figure(1);
subplot(2, 1, 1)
h = imagesc(PETHts, 1:PETHevnr, (PETHdata));
set(h, 'AlphaData', ~isnan(PETHdata))
clim([-5 5])
subplot(2, 1, 2)
g = plot(PETHts, median(PETHdata, 1, 'omitnan'), 'LineWidth', 2);
ylim([-3 3])
xlim([min(PETHts) max(PETHts)])
xline(0, '-r')
yline(mean(mean(BLdataZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)
titlefig1 = string(strcat(combiname,' Event Onset'));
set(fig1, 'Name', titlefig1)

fig2 = figure(2);
plot(Time, contvar)
hold on
for ii = 1:size(Events, 1)
    xline(Events(ii, 1), '-g')
    xline(Events(ii, 2), '-r')
end
title('ZScore DFF with onset freezing events')
titlefig2 = string(strcat(combiname,' ZScore with onset events'));
set(fig2, 'Name', titlefig2)

% Sort freezing episodes by length: 

[~, idx] = sort(boutdur, 'ascend');
%in case of repeated durations:

for  ii= 1:size(boutdur, 1)
    sortedPETH(ii, :) = PETHev((idx(ii)), :);
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
yline(mean(mean(BLdataZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)

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
yline(mean(mean(BLdataZ, 1, 'omitnan'), 'omitnan'), '-', 'LineWidth', 1)


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

PETHfsev = bin;
PETHtsAUC = PETHts + Prewfr; % Escalo a 0

for jj = 1:size(binnedBL, 1)
    if double(sum(~isnan(binnedBL(jj, :))))*PETHfsev < 0.5 % criterio de duración.. si ponemos una fs para el PETH muy alto, hay que tener cuidado porque podría pasar mucho esto al ser un ds
        tbl(jj) = sum(~isnan(binnedBL(jj, :)))*PETHfsev; AUCbl(jj) = nan;
        densbl(jj) = nan; AUCev(jj) = nan; densev(jj) = nan;
    else
    tbl(jj) = double(sum(~isnan(binnedBL(jj, :))))*PETHfsev; % duración bin en seg // la dur del bout ya la tenemos
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
    for jj = 1:size(BLdataZ, 1)
        if double(sum(PETHts >= AUCint(ii, 1) & PETHts < AUCint(ii, 2) & ~isnan(PETHev(jj, :))))*bin < 0.5
            AUCdata(jj, ii) = nan;
            DENSdata(jj, ii) = nan;
            AUCdur(jj, ii) = double(sum(PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHev(jj, :))))*bin;
        else
        AUCdata(jj, ii) = trapz(bin, PETHev(jj, (PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHev(jj, :)))));
        AUCdur(jj, ii) = double(sum(PETHts >= AUCint(ii, 1) & PETHts <= AUCint(ii, 2) & ~isnan(PETHev(jj, :))))*bin;
        DENSdata(jj, ii) = AUCdata(jj, ii)./length(PETHev(jj, (PETHts >= AUCint(ii, 1) & PETHts < AUCint(ii, 2) & ~isnan(PETHev(jj, :)))));
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

prepostFP = table(AUCev.',AUCbl.', densev.', densbl.', tbl.', boutdur, 'VariableNames', measurenames1);


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

%% Calculate n for each value in mean:

idxPETH = ~isnan(PETHev); nint = sum(idxPETH, 1);
idxudPETH = ~isnan(udPETHdata); nintud = sum(idxudPETH, 1);

%% 12. Save all data and figures: 

% Params:
results.PETH.offset.combiname.params.Pre = Pre;
results.PETH.offset.combiname.params.Post = Post;
results.PETH.offset.combiname.params.Prewfr = Prewfr;
results.PETH.offset.combiname.params.PreWind = PreWind;
results.PETH.offset.combiname.params.PostWind = PostWind;
results.PETH.offset.combiname.params.bin = bin;
results.PETH.offset.combiname.params.binsize = binsize; % in samples
results.PETH.offset.combiname.params.binthreshold = threshold;
results.PETH.offset.combiname.params.summarymeasures = summarymeasures;

% PETH variables: 

results.PETH.offset.combiname.PETHdata = PETHfreez;
results.PETH.offset.combiname.PETHwbincriteria = PETHdata;
results.PETH.offset.combiname.PETHts = PETHts;
results.PETH.offset.combiname.udPETHdata = udPETHdata;
results.PETH.offset.combiname.udPETHts = udPETHts;

results.PETH.offset.combiname.PETHmean = PETHm;
results.PETH.offset.combiname.PETHdev = PETHd;
results.PETH.offset.combiname.udPETHmean = udPETHm;
results.PETH.offset.combiname.udPETHdev = udPETHd;

results.PETH.offset.combiname.sortedPETH = sortedPETH;
results.PETH.offset.combiname.EvDFFZ = EvDFFZ;
results.PETH.offset.combiname.BLDFFZ = BLDFFZ;



% Figures: 
% savepath = uigetdir();
folderpath1 = strcat(savepath, '/Figures');
folderpath2 = strcat(savepath, '/Data');
mkdir(folderpath1, 'offsetPETH');
mkdir(folderpath2, 'offsetPETH')

fig1name = strcat(savepath, '/Figures/offsetPETH/offcompletePETH.jpeg');
fig2name = strcat(savepath, '/Figures/offsetPETH/offeventsDFF.jpeg');
fig3name = strcat(savepath, '/Figures/offsetPETH/offsortedPETH.jpeg');
fig4name = strcat(savepath, '/Figures/offsetPETH/offblthreshPETH.jpeg');

saveas(fig1, fig1name)
saveas(fig2, fig2name)
saveas(fig3, fig3name)
saveas(fig4, fig4name)

% Tables:

xlsx1 = strcat(savepath, '/Data/offsetPETH/offprepostAUC.xlsx');
xlsx2 = strcat(savepath, '/Data/offsetPETH/offallintAUC.xlsx');
xlsx3 = strcat(savepath, '/Data/offsetPETH/offsamelengthintAUC.xlsx');
 if exist(xlsx1, 'file')
 delete(xlsx1);
 end
 if exist(xlsx2, 'file')
 delete(xlsx2);
 end
 if exist(xlsx3, 'file')
 delete(xlsx3);
 end
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
xlsx4 = strcat(savepath, '/Data/offsetPETH/completePETH.xlsx');
xlsx5 = strcat(savepath, '/Data/offsetPETH/sortedPETHwthresh.xlsx');
xlsx6 = strcat(savepath, '/Data/offsetPETH/PETHwthresh.xlsx');
PETHev2 = array2table(PETHev, 'VariableNames', pethcolcompletets);
sortedPETH2 = array2table(sortedPETH, 'VariableNames', pethcolcompletets);
PETHdata2 = array2table(PETHdata, 'VariableNames',pethcolcompletets);
 if exist(xlsx4, 'file')
 delete(xlsx4);
 end
 if exist(xlsx5, 'file')
 delete(xlsx5);
 end
 if exist(xlsx6, 'file')
 delete(xlsx6);
 end
writetable(PETHev2, xlsx4, 'WriteVariableNames', true);
writetable(sortedPETH2, xlsx5, 'WriteVariableNames', true);
writetable(PETHdata2, xlsx6, 'WriteVariableNames', true);

pethcolintts = string(udPETHts);
xlsx7 = strcat(savepath, '/Data/offsetPETH/udPETH.xlsx');
udPETHdata2 = array2table(udPETHdata, 'VariableNames', pethcolintts);
 if exist(xlsx7, 'file')
 delete(xlsx7);
 end
writetable(udPETHdata2, xlsx7, 'WriteVariableNames', true);



% Mean + sem:

    
if size(PETHdata, 1) == 1
    PETHd = nan(length(PETHts),1);
    udPETHd = nan(length(udPETHts),1);
    summaryPETH = table(PETHts.', PETHm.', PETHd, nint.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});
    udsummaryPETH = table(udPETHts.', udPETHm.', udPETHd, nintud.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});
else
    summaryPETH = table(PETHts.', PETHm.', PETHd.', nint.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});
    udsummaryPETH = table(udPETHts.', udPETHm.', udPETHd.', nintud.', 'VariableNames', {'Time' char(summarymeasures(1)) char(summarymeasures(2)) 'n'});
end


xlsx8 = strcat(results.FP.path, '/Data/offsetPETH/offsummaryPETH.xlsx');
xlsx9 = strcat(results.FP.path, '/Data/offsetPETH/offudsummaryPETH.xlsx');
if exist(xlsx8, 'file')
delete(xlsx8);
end
writetable(summaryPETH, xlsx8)

if exist(xlsx9, 'file')
delete(xlsx9);
end
writetable(udsummaryPETH, xlsx9)




end