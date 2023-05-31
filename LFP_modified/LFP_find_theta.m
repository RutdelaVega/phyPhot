function [ theta_analysis] = LFP_find_theta(basicData,varargin)
% Perform theta analysis. Compute phase, raleight and mvl
%
%   INPUTS
%       basicData (struct): Structure containing basic neural data
%       config_channels (struct): Structure containing the LCN_set_config_channels parameters
%
%   OPTIONAL INPUT
%       fsLow (int): to download (4800 by default)
%       spec (array): spectral band ([4 12] for theta by default)
%       idxsDiscard (array): indexes of lfp not to take to analysis
%       peak_trough (char): explicit if peak to peak, or trough to trough (by default that 
%                          so in pyramidal is trough to trough)
%       consecutive (logic): if true, erase all cycles that are not surrounded by other theta cycles
%       filter_type (char): explicit type of filter (by default butter)
%       do_cyclemat (logic): if true, make matrix with theta cycles
%
%   OUTPUTS
%        theta_analysis (struct): Structure containing theta analysis
%
%     Some info about the output
%           .iCyc: indexes of beginning and end of cycles
%           .iMiddle: indexes of middle of cycles
%           .spkPhases: phases of spikes
%           .ph_spec: structure with output data:
%               .phase: mean phase (in deg)
%               .sd: phase standard deviation (in deg)
%               .mvl: mean vector lenght
%               .Z: Z Rayleigh's test
%               .p: p Rayleigh's test
%           .dataHil: hilbert of filtered data
%           .dataFilt: filtered data
%
% LCN - MV 2018
% Modified: A Navas-Olive 2019
% Modified: Cecilia Gallego-Carracedo 2022


    % Parse inputs 
    p = inputParser;
    addParameter(p,'fsLow',4800,@isnumeric);
    addParameter(p,'spec',[4 12],@ismatrix);
    addParameter(p,'idxsDiscard',[],@ismatrix);
    addParameter(p,'peak_trough','',@ischar);
    addParameter(p,'consecutive',0,@islogical);
    addParameter(p,'filter_type','butter',@ischar);
    addParameter(p,'do_cyclemat', true, @islogical);
    addParameter(p,'min_undiscard', 2, @isnumeric); % seconds
    parse(p,varargin{:});
    fsLow = p.Results.fsLow;
    spec = p.Results.spec;
    idxsDiscard = p.Results.idxsDiscard;
    peak_trough = p.Results.peak_trough;
    consecutive = p.Results.consecutive;
    filter_type = p.Results.filter_type;
    do_cyclemat = p.Results.do_cyclemat;
    min_undiscard = p.Results.min_undiscard;

    LFP = basicData.dataLFP; 
    fs = basicData.fsLFP; 

%    areaTheta = basicData.area;
    
    fprintf('\n\n    --------------------\n    |  Theta Analysis  |\n    --------------------\n');
    

    % Discard times
    undiscarded_seconds = sum(~idxsDiscard)/fs;
    if undiscarded_seconds < min_undiscard
        idxsDiscard = false(size(LFP,1),1);
    end    

    % Downsample    
    iLong = 1 : size(LFP,1);                                % Indexes for long LFP
        fsLow = fs;
        iShort = iLong;
        data = LFP;  
        idxs2nan = round(find(idxsDiscard));
    idxs2nan(idxs2nan<1) = 1; idxs2nan(idxs2nan>size(data,1)) = size(data,1);
    
    % Peak or trough?

    % Automatic search for taking trough to trough in SP

    % Specified by user
    if strcmp(peak_trough, 'peak')
            PT = 'peak';
            disp('     > trough to trough in SP, by peak to peak in SR/SLM')
    elseif strcmp(peak_trough, 'trough')
            PT = 'trough';
            disp('     > trough to trough in SP, by trough to trough in SP/SO')
    end


    % Filter signal
    switch filter_type

        case 'firls'

            % Filter
            N      = 10000;  % Order, works as a windowing
            Fstop1 = spec(1)-.1;   % First Stopband 	
            Fpass1 = spec(1);   % First Passband Frequency
            Fpass2 = spec(2);  % Second Passband Frequency
            Fstop2 = spec(2)+.1;  % Second Stopband Frequency
            Wstop1 = 1;     % First Stopband Weight
            Wpass  = 1;     % Passband Weight
            Wstop2 = 1;     % Second Stopband Weight
            % Filtering
            bRP  = firls(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 fsLow/2]/(fsLow/2), [0 0 1 1 0 0], [Wstop1 Wpass Wstop2]);
            dataFilt = filtfilt(bRP,1,double(data));

        case 'butter' 

            N = 2;
            Fpass1 = spec(1);
            Fpass2 = spec(2);

            [b, a] = butter(N, [Fpass1 Fpass2]/(fs/2));
            dataFilt = filtfilt(b, a, double(data)); 

    end



    % Compute phases from dataFilter (before it has nans)
    switch PT
        case 'peak'
            % Hilbert: analitical signal (in complex space)
            dataHil = hilbert(dataFilt); 
        case 'trough'
            % Hilbert: analitical signal (in complex space)
            dataHil = hilbert(-dataFilt); 
    end
        
    % Erase parts of LFP
    if ~isempty(idxs2nan)
        dataFilt(idxs2nan) = nan;
    end

    switch PT
        case 'peak' % Find peaks
            [ ~, locs ] = findpeaks( + dataFilt ); 
        case 'trough' % Find troughs
            [ ~, locs ] = findpeaks( - dataFilt );
    end

    % Make matrix of those cycle edges
    iCyc = [ locs(1:end-1)';... % Beginning of cycles
             locs(2:end)' ];  % End of cycles
    %iCyc = iCyc./fs;

    % Length of cycles
    Lcyc = iCyc(2,:) - iCyc(1,:);
    % Keep only those that are in the spectral band:
    iCyc = iCyc( :, Lcyc./fs >= 1/spec(2) & Lcyc./fs <= 1/spec(1) );
    % Find peak/troughs
    iMiddle = zeros(1,size(iCyc,2));
    goodMiddles = ones(1,size(iCyc,2));
    for ii = 1:length(iMiddle)
        idxs = iCyc(1,ii):iCyc(2,ii);
        % Find peaks
        switch PT
            case 'peak'
                [ ~, imax ] = nanmin( dataFilt( idxs ) );
            case 'trough'
                [ ~, imax ] = nanmax( dataFilt( idxs ) ); 
        end
        imid = idxs(1) + imax -1;
        if imid >= mean(idxs) - 0.3*std(idxs) && imid <= mean(idxs) + 0.3*std(idxs)
            iMiddle(ii) = imid;
        else
            goodMiddles(ii) = 0;
        end
    end
    % Discard bad middles, and those cycles
    iCyc = iCyc( :, goodMiddles==1 );
    iMiddle = iMiddle( goodMiddles==1 );
    
    % Hights of theta cycles
    hights = [];
    for ii = 1:length(iMiddle)
        hights =[ hights, max(data(iCyc(1,ii):iCyc(2,ii)))-min(data(iCyc(1,ii):iCyc(2,ii))) ];
    end
    
    figure('pos',[100,100,800,500]); hold on
    [ yhist, xhist ] = hist( double(hights), 50 );
    bar(xhist,yhist,1,'facecolor',[.6 .6 .6])
    xmax = xhist(yhist==max(yhist)); xmax = xmax(1);
    dragbarLow = plot( xmax-[1 1], [0 max(yhist)*1.1], 'k', 'linewidth', 3 );
    draggable(dragbarLow,'constraint','h');
    dragbarHigh = plot( xhist(end)+[1 1], [0 max(yhist)*1.1], 'k', 'linewidth', 3 );
    draggable(dragbarHigh,'constraint','h');
    xlim([xhist(1)-2,xhist(end)+2])
    xlabel('Theta cycles hights')
    title('Select minimum theta cycles hights')
    btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
        'Units','normalize','Position', [.73 .72 .09 .06],...
        'Callback', 'uiresume(gcbf)');
    %uiwait(gcf);
    get(dragbarLow,'XData');
    minHight = mean(get(dragbarLow,'XData'));
    get(dragbarHigh,'XData');
    maxHight = mean(get(dragbarHigh,'XData'));
    %close;
    
    % Discard cycles over maxHight or lower than minHight
    iCyc = iCyc( :, (hights>minHight) & (hights<maxHight) );
    iMiddle = iMiddle( :, (hights>minHight) & (hights<maxHight) );
    
    % Discard not-consecutive cycles
    if consecutive
        time_around = [iCyc(1,2:end) - iCyc(2,1:end-1), 0];
        consecutive_cycles = 1+find( (time_around(1:end-1)/fs <= 0.2) & (time_around(2:end)/fs <= 0.2) );
        iCyc = iCyc(:, consecutive_cycles);
        iMiddle = iMiddle(consecutive_cycles);
    end

    % Upsample the output data
    iCyc = round(iCyc / fsLow*fs);                     % Indexes where theta cycles start (first row), and end (second row)
    iMiddle = round(iMiddle / fsLow*fs);               % Indexes where theta maximum is
    dataFilt = interp1(iShort, dataFilt, iLong )';    % Upsample LFP

    
    % Outputs
    theta_analysis.iCyc = iCyc;
    theta_analysis.iMiddle = iMiddle;
    theta_analysis.chTheta = LFP;

    % Optional outputs
    theta_analysis.dataHil = dataHil;
    theta_analysis.dataFilt = dataFilt;


end