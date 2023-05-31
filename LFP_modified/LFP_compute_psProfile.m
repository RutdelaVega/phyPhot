function [psProfile] = LFP_compute_psProfile(basicData, varargin)
% Power channel profile and power spectrum over windows + hanning
%   INPUTS
%       basicData (struct): Structure containing basic neural data in different entries
%       config_channels (struct): Structure containing the LCN_set_config_channels parameters
%
%   OPTIONAL INPUTS
%       saveName (char): Name to save output image
%       tWin (numeric): Window size (seconds) to compute the spectrogram
%                       (default: 2s)
%       dtWin (numeric): step size (seconds) of each spectrogram window
%                       (default: 0.5s)        
%   OUTPUTS
%       pS (matrix) 3d power spectrum matrix (frequencies x time windows (2s) x channel)
%
% LCN-MV 2016
% Edit: LCN-AN 2019
% Edit: Cecilia Gallego-Carracedo 2022
% Edit: Julio Esparza 13-05-2022

% Parse inputs 
p = inputParser;
addParameter(p,'saveName','',@isstr); %#ok<DISSTR> 
addParameter(p, 'tWin', 2, @isnumeric); % 2 second window
addParameter(p, 'dtWin', 0.5, @isnumeric); % 0.5 seconds side step

parse(p,varargin{:});
saveName = p.Results.saveName;
tWin = p.Results.tWin;
dtWin = p.Results.dtWin;

data = basicData.dataLFP;
fs = basicData.fsLFP;
%shMap = config_channels.shMap; 
%shType = config_channels.shType;

delta_cf = [0, 3]; % Delta power: between 0 and 3Hz
theta_cf = [4, 12]; % Theta power: between 4 and 12Hz
alpha_cf = [9, 16]; % Alpha power: between 9 and 16Hz
gamma_cf = [20, 90]; % Gamma power: between 20 and 90Hz
hfo_cf = [100, 1e3]; % HFO power: between 100 and 1000Hz
ch = 1;

fprintf('\n\n    --------------------\n    |  Power Spectrum  |\n    --------------------\n');
    % Electrode position, mm with respect bregma

    % Data matrix needs to be time x channels
    if size(data,2) > size(data,1) 
        data = data';
    end

    keepInd = 1:size(data,2);
    data = data - cast(mean(data,1, 'omitnan'), 'like', data);

    % DO POWER SPECTRUM IN WINDOWS OF 2s EVERY 0.5s FOR EACH CHANNEL
    tWin = fs*tWin;
    % Window moving time step
    dtWin = fs*dtWin;
    % Number of windows
    nWin = floor(size(data,1)/dtWin) - 3;
    % Power Spectrum storing matrix
    lenpS = length(pspectrum(double(data(1:round(tWin), 1)'),  fs, 'Reassign', true)); % Por qué coge sólo la primera parte de los datos? Cuando saca el pspectrum se dobla
    pS = zeros(lenpS,nWin,1);

    % For each channel
    %for ch = channels % Recorro todos los canales.
        for win = 1:nWin
            % Temporal data of the corresponding window
            dataWin = data(round(dtWin*(win-1)+1):round(dtWin*(win-1)+tWin),1);
            % Storing results of power spectrum between 1 and 1000
            [pS(:,win,ch),freq] = pspectrum(double(dataWin'), fs, 'Reassign', true);
        end
    %end

    % For each channel
    pdelta = nan(nWin, 1);
    ptheta = nan(nWin, 1);
    palpha = nan(nWin, 1);
    pgamma = nan(nWin, 1);
    phfo   = nan(nWin, 1);

    %for ch = channels

        % For each window
        for win = 1:nWin
            % Delta power: between 0 and 3Hz
            pdelta(win) = 10*log10(sum(pS(freq>=delta_cf(1)&freq<=delta_cf(2), win)))+60;
            % Theta power: between 4 and 12Hz
            ptheta(win) = 10*log10(sum(pS(freq>=theta_cf(1)&freq<=theta_cf(2), win)))+60;
            % Alpha power: between 9 and 16Hz
            palpha(win) = 10*log10(sum(pS(freq>=alpha_cf(1)&freq<=alpha_cf(2), win)))+60;
            % Gamma power: between 20 and 90Hz
            pgamma(win)  = 10*log10(sum(pS(freq>=gamma_cf(1)&freq<=gamma_cf(2), win)))+60;
            % HFO power: between 100 and 1000Hz
            try %y este catch, cuando sucede??
                phfo(win) = 10*log10(sum(pS(freq>=hfo_cf(1)&freq<=hfo_cf(2), win)))+60; 
            catch 
                phfo(win) = nan; 
            end 
        end
    %end
    % Means and stds for each rhythm
    meanDelta = mean(pdelta);       stdDelta = std(pdelta);     chDelta = find(meanDelta ~= -Inf );
    meanTheta = mean(ptheta);       stdTheta = std(ptheta);     chTheta = find(meanTheta ~= -Inf);
    meanAlpha = mean(palpha);       stdAlpha = std(palpha);     chAlpha = find(meanAlpha ~= -Inf );
    meanGamma = mean(pgamma);       stdGamma = std(pgamma);     chGamma = find(meanGamma ~= -Inf);
    meanHfo   = mean(phfo);         stdHfo   = std(phfo);       chHfo   = find(meanHfo ~= -Inf);
    % channels without nans
    chNoNaN = find(any(~isnan([meanHfo; meanGamma; meanTheta; meanDelta; meanAlpha]),1));
    chDelta = chDelta(chNoNaN);
    chTheta = chTheta(chNoNaN);
    chAlpha = chAlpha(chNoNaN);
    chGamma = chGamma(chNoNaN);
    chHfo   = chHfo(chNoNaN);

    % MAKE FIGURE
    fig1 = figure('pos',[100,10,40*size(data,2),800]); hold on;
    ymax = 1.3*max( [ meanHfo meanGamma meanTheta meanDelta meanAlpha ] );
    ymin = 0.8*min( [ meanHfo meanGamma meanTheta meanDelta meanAlpha ] );
    plot( keepInd([ 1:size(data,2) ; 1:size(data,2) ]), [ymin*ones(1,size(data,2)); ymax*ones(1,size(data,2))], '-', 'color', [.8 .8 .8] )

    % Plot std with shadows
    fill(keepInd([chHfo flip(chHfo) chHfo(1)]),...
        [meanHfo(chHfo)+stdHfo(chHfo) flip(meanHfo(chHfo)-stdHfo(chHfo))...
        meanHfo(chHfo(1))+stdHfo(chHfo(1))],[.1 .1 .1],'facealpha',.1,'lineStyle','none');
    fill(keepInd([chGamma flip(chGamma) chGamma(1)]),...
        [meanGamma(chGamma)+stdGamma(chGamma) flip(meanGamma(chGamma)-stdGamma(chGamma))...
        meanGamma(chGamma(1))+stdGamma(chGamma(1))],[.2 .9 .2],'facealpha',.1,'lineStyle','none');
    fill(keepInd([chTheta flip(chTheta) (chTheta(1))]),...
        [meanTheta(chTheta)+stdTheta(chTheta) flip(meanTheta(chTheta)-stdTheta(chTheta))...
        meanTheta(chTheta(1))+stdTheta(chTheta(1))],[.7 .2 .2],'facealpha',.1,'lineStyle','none');
    fill(keepInd([chDelta flip(chDelta) (chDelta(1))]),...
        [meanDelta(chDelta)+stdDelta(chDelta) flip(meanDelta(chDelta)-stdDelta(chDelta))...
        meanDelta(chDelta(1))+stdDelta(chDelta(1))],[.2 .2 .7],'facealpha',.1,'lineStyle','none');
    fill(keepInd([chAlpha flip(chAlpha) (chAlpha(1))]),...
        [meanAlpha(chAlpha)+stdAlpha(chAlpha) flip(meanAlpha(chAlpha)-stdAlpha(chAlpha))...
        meanAlpha(chAlpha(1))+stdAlpha(chAlpha(1))],[.7 .2 .7],'facealpha',.1,'lineStyle','none');
    % Plot means
    i1 = plot( keepInd(chHfo), meanHfo(chHfo), '-o', 'color', [.1 .1 .1], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
    i2 = plot( keepInd(chGamma), meanGamma(chGamma), '-o', 'color', [.2 .9 .2], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
    i3 = plot( keepInd(chTheta), meanTheta(chTheta), '-o', 'color', [.7 .2 .2], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
    i4 = plot( keepInd(chDelta), meanDelta(chDelta), '-o', 'color', [.2 .2 .7], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );
    i5 = plot( keepInd(chAlpha), meanAlpha(chAlpha), '-o', 'color', [.7 .2 .7], 'LineWidth' , 2, 'Markersize', 8, 'MarkerFaceColor', 'w' );

    % Legend
    ylabel('Power (dB)'); xlabel('# Channel');
    legend([i1 i2 i3 i4 i5],'HFOs 100-1000 Hz','Gamma 20-90 Hz','Theta 4-12 Hz',...
        'Alpha 9-15 Hz','Delta 1-3 Hz','location','southeast');
    xlim([0 keepInd(size(data,2))+.5]);
    ylim([min([meanHfo(chHfo) meanGamma(chGamma) meanTheta(chTheta) meanDelta(chDelta) meanAlpha(chAlpha)])*0.99, max([meanHfo(chHfo) meanGamma(chGamma) meanTheta(chTheta) meanDelta(chDelta) meanAlpha(chAlpha)])*1.05])
    xTickVal = keepInd(1:2:size(keepInd, 2));
    set(gca,'xtick',xTickVal)
    % Separate shanks
    ax = axis;
    plot( [.5 .5], [ax(3) ax(4)*1.02], '--', 'color', [.5 .5 .5], 'linewidth', 2 )
    hold off


    % IF timeSpec IS TRUE, PLOTS THE POWER SPECTRUM OF THETA
    % Question: do you want to see the power Spectrum alon time for the channel that has more theta power?
    qu = 'n';
    % If you think is correct, then 'yes', to distinguish between theta and LIA
    if strcmp('y',qu) || strcmp('1',qu)
        % Computes channel with maximum Theta
        qu2 = input(sprintf('Channel with max theta is %d, do power spectrum (y/n)?\n',find(max(meanTheta)==meanTheta)),'s');
        if strcmp('y',qu2) || strcmp('1',qu2)
            chSLM = find(max(meanTheta)==meanTheta);
        else
            chSLM = input('Channel for power spectrum?: ');
        end
    else
        chSLM = nan;
    end

    if ~isnan(chSLM)
        % En cada tags calculo el power de theta (columna 3 de tags) y de delta (columna 4 de tags)
        % y hago el cociente (Senior et al., 2008) (columna 5 de tags)
        % Power spectrum of SLM in decibels
        pSslm = 10*log10(pS(:,:,chSLM))+60;
        % Compute relation between theta and LIA powers (Senior et al., 2008)
        tvsd = nan(nWin,3);
        for win = 1:nWin
            tvsd(win,1) = sum(pSslm(freq>=5 & freq<=9, win));
            tvsd(win,2) = sum(pSslm(freq>=1 & freq<=3, win));
            tvsd(win,3) = tvsd(win,1)/tvsd(win,2);
        end

        fig2 = figure; % Hago un espectograma y le ploteo encima la relaci�n theta/delta
        subplot(7,1,2:7)
        imagesc([1 length(data)/fs],(freq(freq>=1 & freq<=30)),(pSslm((find(freq>=1 & freq<=30)),:)));
        xlabel('s)');
        ylabel('freq(Hz)');
        colormap jet
        h11 = colorbar('southoutside'); h11.Label.String = 'Power (dB)';
        subplot(7,1,1)
        plot(tvsd(:,3));
        xlim([1 length(tvsd)]); ylim([1 3]);
        ylabel('th/del');
    end
    
    % Plot
    if ~isempty(saveName)
        saveas(gcf, saveName);
    end

    % Saving info
    % Theta
    band.theta.pmean = meanTheta;
    band.theta.pstd = stdTheta;
    band.theta.ch = chTheta;
    % Delta
    band.delta.pmean = meanDelta;
    band.delta.pstd = stdDelta;
    band.delta.ch = chDelta;
    % Gamma
    band.gamma.pmean = meanGamma;
    band.gamma.pstd = stdGamma;
    band.gamma.ch = chGamma;
    % HFO
    band.hfo.pmean = meanHfo;
    band.hfo.pstd = stdHfo;
    band.hfo.ch = chHfo;
    % Alpha
    band.alpha.pmean = meanAlpha;
    band.alpha.pstd = stdAlpha;
    band.alpha.ch = chAlpha;
    % Params
    params = struct();
    params.delta_cf = delta_cf;
    params.theta_cf = theta_cf;
    params.alpha_cf = alpha_cf;
    params.gamma_cf = gamma_cf;
    params.hfo_cf = hfo_cf;
    params.tWin = tWin;
    params.dtWin = dtWin;

    % Save all in meta-structure to be outputed
    psProfile.pS = pS;
    psProfile.infoRhs = band;
    psProfile.pSfreqs = freq;
    psProfile.handle_Profile = fig1;
    psProfile.params = params;
    if exist('fig2','var'); psProfile.handle_pS = fig2; end
end
