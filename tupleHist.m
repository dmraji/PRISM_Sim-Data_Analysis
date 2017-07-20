function [smEnTC, timeAdj] = tupleHist()
    
    % NOTE: File read-in currently does not work with Octave

    % Constants for decay (Cs-137)
    % Assuming source strength of 10 microCi
    lam = 0.693 / (30 * 3.154 * (10^7));
    N = (10.08 * 3.7*(10^4)) / lam;
    act = lam * N / 1000000000;

    % Grabbing data from G4 PRISM_Sim output file, transforming into matrix
    % pr = 'Enter file name.';
    % fn = input(pr, 's');
    fn = 'output_662keV_1Det_coneWide_phi0_theta90.txt';
    fID = fopen(fn, 'r');
    line = fgetl(fID);
    fclose(fID);
    isStrCol = isnan(str2double(regexp(line, '[^\t]+', 'match')));
    
    format = cell(1, numel(isStrCol));
    format(isStrCol) = {'%s'};
    format(~isStrCol) = {'%f'};
    format = [format{:}];
    
    fID = fopen(fn, 'r');
    data = textscan(fID, format, Inf, 'Delimiter', '\t');
    fclose(fID);
    
    for colId = find(~isStrCol)
        data{colId} = num2cell(data{colId});
    end
    data = [data{:}];
    dataNC = cellfun(@str2double, data);
    
    % Grabbing individual vectors from data matrix
    evNum = dataNC(2:end, 1);
    hitNum = dataNC(2:end, 2);
    trackID = dataNC(2:end, 3);
    en = dataNC(2:end, 4);
    detID = dataNC(2:end, 5);
    process = dataNC(2:end, 6);
    DOI = dataNC(2:end, 7);
    Phi = dataNC(2:end, 8);
    Theta = dataNC(2:end, 9);
    HP = dataNC(2:end, 10);
    time = dataNC(2:end, 11);
    
    % Adjustments for "noise" from sim - largely result of ...
    % transport "scatterings" in G4
    noNoiseEn = en;
%     histogram(noNoiseEn, 256)
%     title('Energy, Binned');
%     xlabel('Energy, keV');
%     ylabel('Counts');
%     grid on;
%     fprintf('Press any key to continue.\n');
%     pause
    
%     histogram(DOI, 10)
%     title('Depth of Interaction, Binned');
%     xlabel('Depth of Interaction (mm)');
%     ylabel('Counts');
%     grid on;
    
    % Coincidence considerations
    i = 1;
    ii = 1;
    
    coinTr = evNum;
    timeAdj = time;
    timeSaver = 0;
    TCen = en;
    
    % Considering events that miss the detector geometry, adjusting ...
    % timing accordingly
    while ii <= evNum(end) && i <= length(evNum)
        if ~any(i==evNum)
            while ~any(ii==evNum)
                 timeSaver = timeAdjGhost(timeSaver);
                 ii = ii + 1;
            end
        end
        if i == 1
            coinTr(i) = 0;
            timeSaver = timeAdjUq(i, timeSaver);
            i = i + 1;
            ctr = 1;
            while evNum(i-1) == evNum(i)
                tTracker = time(i) - time(i-1);
                
                coinTr(i) = evNum(i);
                coinTr(i-1) = evNum(i-1);
                
%                 TCen(i-ctr) = TCen(i-ctr) + TCen(i);
%                 TCen(i) = 0;
                timeAdj(i) = timeAdj(i-1) + tTracker;
                i = i + 1;
                ctr = ctr + 1;
            end
        else
            ctr = 1;
            try
                while evNum(i-1) == evNum(i)
                    tTracker = time(i) - time(i-1);

                    coinTr(i) = evNum(i);
                    coinTr(i-1) = evNum(i-1);

%                     TCen(i-ctr) = TCen(i-ctr) + TCen(i);
%                     TCen(i) = 0;

                    timeAdj(i) = timeAdj(i-1) + tTracker;
                    i = i + 1
                    ctr = ctr + 1;
                end
            end
            if evNum(i-1) ~= evNum(i)
                coinTr(i) = 0;
                timeSaver = timeAdjUq(i, timeSaver);
                i = i + 1;
                try
                    ctr = 1;
                    while evNum(i-1) == evNum(i)
                        tTracker = time(i) - time(i-1);
                        
                        coinTr(i) = evNum(i);
                        coinTr(i-1) = evNum(i-1);
                        
%                         TCen(i-ctr) = TCen(i-ctr) + TCen(i);
%                         TCen(i) = 0;
                        
                        timeAdj(i) = timeAdj(i-1) + tTracker;
                        i = i + 1;
                        ctr = ctr + 1;
                    end
                end
            end
        end
        ii = ii + 1;
    end
    
    coinEn = coinTr .* en;
    coinEn = coinEn(find(coinTr));
    
    % Timing adjustment for lack of true "global" time in simulation
    % Also adjusting energy vector due to limitations of real detector ...
    % (summing multi-track events contained within detector)
    function [timeSaver] = timeAdjUq(ind, timeSaver)
        timePl = -1 * log(rand) / (act);
        if timeSaver == 0
            timeSaver = timePl;
        else
            timePl = timePl + timeSaver;
            timeSaver = timePl;
        end
        timeAdj(ind) = timeAdj(ind) + timePl;
    end

    function [timeSaver] = timeAdjGhost(timeSaver)
        timePl = -1 * log(rand) / (act);
        if timeSaver == 0
            timeSaver = timePl;
        else
            timePl = timePl + timeSaver;
            timeSaver = timePl;
        end
    end

    % identifying "false" coincidences based on time (including ...
    % adjustments)
    % Also adjusting energy vector to sum chance coincident depositions ...
    % and for charge loss with DOI consideration
    
    muTaoEl = 0.01;
    muTaoHo = 0.0008;

    bias = 1000;

    detdepth = 1;
    
    k = 1;
    coinAll = evNum;
    while k < length(timeAdj)
        % Note: 1.3 mircoseconds discussed in meeting as being current
        % goal = 0.45 microSec
        if ((timeAdj(k) + 1000) >= timeAdj(k+1)) && ((timeAdj(k) - 1000) <= timeAdj(k+1))
            
            coinAll(k) = 1;
            coinAll(k+1) = 1;
            
            TCen(k) = TCen(k) + TCen(k+1);
            TCen(k+1) = 0;
            DOI(k) = -1;
            doiHold = DOI(k+1);
            DOI(k) = doiHold;
            DOI(k+1) = -1;
            ctr1 = 2;
            try
                while ((timeAdj(k) + 1000) >= timeAdj(k+ctr1)) && ((timeAdj(k) - 1000) <= timeAdj(k+ctr1))
                    coinAll(k+ctr1) = 1;
                    TCen(k) = TCen(k) + TCen(k+ctr1);
                    TCen(k+ctr1) = 0;
                    DOI(k+ctr1-1) = -1;
                    doiHold = DOI(k+ctr1);
                    DOI(k) = doiHold;
                    DOI(k+ctr1) = -1;
                    ctr1 = ctr1 + 1;
                end
            end
            
            % Charge collection efficiency
%             CCE = ((muTaoEl * bias)/(detdepth) * (1 - exp(-(detdepth - (DOI(k) * 0.1)) / (muTaoEl * bias)))) + ((muTaoHo * bias)/(detdepth) * (1 - exp(-(DOI(k) * 0.1) / (muTaoHo * bias))));
%                 
%             enCl = TCen(k) .* CCE;
%             TCen(k) = enCl;
            
            k = k + ctr1;
        else
            coinAll(k) = 0;
            
%             CCE = ((muTaoEl * bias)/(detdepth) * (1 - exp(-(detdepth - (DOI(k) * 0.1)) / (muTaoEl * bias)))) + ((muTaoHo * bias)/(detdepth) * (1 - exp(-(DOI(k) * 0.1) / (muTaoHo * bias))));
%                 
%             enCl = TCen(k) .* CCE;
%             TCen(k) = enCl;
            
            k = k + 1;
        end
    end
    if k == length(timeAdj)
        if en(k) == TCen(k)
            
%             CCE = ((muTaoEl * bias)/(detdepth) * (1 - exp(-(detdepth - (DOI(k) * 0.1)) / (muTaoEl * bias)))) + ((muTaoHo * bias)/(detdepth) * (1 - exp(-(DOI(k) * 0.1) / (muTaoHo * bias))));
%                 
%             enCl = TCen(k) .* CCE;
%             TCen(k) = enCl;
            
        end
    end
    if coinAll(end) ~= 1
        coinAll(end) = 0;
    end
    
    coinAllEn = coinAll .* en;
    coinAllEn = coinAllEn(find(coinAll));
    
    histogram(timeAdj, 100)
    title('Time, Binned');
    xlabel('Time, ns');
    ylabel('Counts');
    grid on;
    fprintf('Press any key to continue.\n');
    pause
    
    % Adjusting data vectors to remove transport processes
    noTrTCen = nonzeros(TCen);
%     noTrDOI = DOI(DOI >= 0);
    
    helpVec = TCen;
    itn = 1;
    while itn <= length(helpVec)
        if helpVec(itn) ~= 0
            helpVec(itn) = 1;
        end
        itn = itn + 1;
    end
    noTrevNum = evNum(find(helpVec));
    noTrhitNum = hitNum(find(helpVec));
    noTrtrackID = trackID(find(helpVec));
    noTrdetID = detID(find(helpVec));
    noTrprocess = process(find(helpVec));
    noTrDOI = DOI(find(helpVec));
    noTrPhi = Phi(find(helpVec));
    noTrTheta = Theta(find(helpVec));
    noTrHP = HP(find(helpVec));
    noTrtime = time(find(helpVec));
    
    tabbedEn = tabulate(round(noTrTCen));

    [noiseFWHM, howTight] = ROITC;
    
    function [noiseFWHM, howTight] = ROITC()
        [peakCts, peakInd] = max(tabbedEn(:, 2));
        peakEn = tabbedEn(peakInd, 1);
        n = 1;
        
        % Assuming a measurement of 9.5 keV-FWHM for electronic noise
        elecNoiseFWHM = 9.5;
        elecNoiseSTDDEV = elecNoiseFWHM / 2.355;
        
        % Statistical noise for peak
        qCarrGenN = peakEn / 0.005;
        fano = 0.1;
        statResLim = 2.355 * sqrt(fano / qCarrGenN);
        statFWHM = statResLim * peakEn;
        statNoiseSTDDEV = statFWHM / 2.355;
        
        % Adding noise stddevs in quadrature
        noiseDev = sqrt((elecNoiseSTDDEV ^ 2) + (statNoiseSTDDEV ^ 2));
        noiseFWHM = 2.355 * noiseDev;
        
        while n <= peakCts
            vecROI(n) = peakEn;
            n = n + 1;
        end
        nFilter = [];
        tightness = 128;
        
        [nFilter, smVecROI, howTight] = tightenerTC(nFilter, tightness);
        
        % Adjusting the fwhm of the spectrum to match that of the ...
        % target (the noise fwhm)
        function [nFilter, smVecROI, howTight] = tightenerTC(nFilter, howTight)
            nFilter = [];
            while length(nFilter) ~= length(vecROI)
                
                filterFiller = randn / howTight; 
                % tightening spread of randn takes resolution at FWHM to  ...
                % approximately 2-percent (13.24 keV is exactly 2%) ...
                % initially, will be revised in recursion later

                if filterFiller
                    nFilter = [nFilter, filterFiller + 1];
                end
            end
            smVecROI = nFilter .* vecROI;
        end
        
        [binInds, binEdges] = discretize(smVecROI, round(length(unique(smVecROI)) / 8));
        tabbedBinInds = tabulate(binInds);
        p = 1;
        halfDiff = 1000;
        diffr = 0;
        fwCur = 1000;
        [rows, cols] = size(tabbedBinInds);
        fact = 1.2;
        upCount = 0;
        tic;
        
        % Tracking down the most accurate bins for the half-max
        while round(noiseFWHM, 1) ~= round(fwCur, 1)
            
            [binInds, binEdges] = discretize(smVecROI, round(length(unique(smVecROI)) / 16));
            tabbedBinInds = tabulate(binInds);
            p = 1;
            halfDiff = 1000000;
            diffr = 0;
            fwCur = 1000000;
            [rows, cols] = size(tabbedBinInds);
            hm = max(tabbedBinInds(:, 2)) / 2;
            
            % two while loops to grab energies of the left and right ...
            % edges of the fwhm
            while p < round(rows / 2)
                diffr = abs(hm - tabbedBinInds(p, 2));
                if diffr < halfDiff
                    halfDiff = diffr;
                    fwIndLeft = (binEdges(p) + binEdges(p + 1)) / 2;
                end
                p = p + 1;
            end

            halfDiff = 1000000;
            while p < rows
                diffr = abs(hm - tabbedBinInds(p, 2));
                if diffr < halfDiff
                    halfDiff = diffr;
                    fwIndRight = (binEdges(p) + binEdges(p + 1)) / 2;
                end
                p = p + 1;
            end

            % Checking whether the current FWHM is smaller or larger ...
            % than that dictated by the noise
            fwCur = fwIndRight - fwIndLeft;
            
            if fwCur > noiseFWHM
                tightness = tightness * fact;
                upCount = upCount + 1;
                if rem(upCount, 5) == 0
                    fact = sqrt(fact);
                end
            elseif fwCur < noiseFWHM
                tightness = tightness / fact;
            end
            
            % Ensuring the recursion doesn't get stuck
            elT = toc;
            if elT > 25
                fact = 1.2;
                tic;
            end
            [nFilter, smVecROI, howTight] = tightenerTC(nFilter, tightness);
        end
        
        histfit(smVecROI, round(length(unique(smVecROI)) / 16))
        title('Fitted peak with randomly-selected normal smearing');
        xlabel('Energy, keV');
        ylabel('Counts');
        grid on;
        fprintf('Press any key to continue.\n');
        pause
        histSmearedTC(peakEn, smVecROI, howTight, elecNoiseSTDDEV);
%         resp = input('Is there another peak? (y / n)\n', 's');
%         if strcmp(resp, 'y')
%             try
%                 j1 = 1;
%                 while j1 <= 9
%                     tabbedEn(peakInd + j1 - 5, 2) = 0;
%                     j1 = j1 + 1;
%                 end
%             end
%             ROITC;
%         elseif strcmp(resp, 'n')
%             fprintf('Moving on ...\n');
%         end
    end
        
    function histSmearedTC(peakEn, smVecROI, howTight, elecNoiseSTDDEV)
        
        x = 1;
        y = 1;
        smEnTC = noTrTCen;
                
        while x <= length(noTrTCen)
            if round(smEnTC(x)) == peakEn
                smEnTC(x) = smVecROI(y);
                y = y + 1;
            else
                
                qCarrGenN = smEnTC(x) / 0.005;
                fano = 0.1;
                statResLim = 2.355 * sqrt(fano / qCarrGenN);
                statFWHM = statResLim * smEnTC(x);
                statNoiseSTDDEV = statFWHM / 2.355;
                
                noiseDev = sqrt((elecNoiseSTDDEV ^ 2) + (statNoiseSTDDEV ^ 2));
                
                noiseFilter = randn * noiseDev;
                smEnTC(x) = smEnTC(x) + (noiseFilter);
            end
            x = x + 1;
            
        end
        
        smEnTC = chargeLoss(smEnTC);
        
        smEnTC = background(smEnTC);

        function [smEnTC] = chargeLoss(smEnTC)
            
            % Numbers based on Michelle's dissertation
            muTaoEl = 0.01;
            muTaoHo = 0.0008;
            
            bias = 1000;
            
            detdepth = 1;
            
            clct = 1;
            while clct <= length(smEnTC)
                % Max CCE is 0.9516, min is 0.5708 (for DOI = 0 & 10, ...
                % resp.)
                CCE = ((muTaoEl * bias)/(detdepth) * (1 - exp(-(detdepth - (noTrDOI(clct) * 0.1)) / (muTaoEl * bias)))) + ((muTaoHo * bias)/(detdepth) * (1 - exp(-(noTrDOI(clct) * 0.1) / (muTaoHo * bias))));
                
                enCL = smEnTC(clct) .* CCE;
                smEnTC(clct) = enCL;
                clct = clct + 1;
            end
            
        end
        
        
        % Adding the measured background into the spectrum
        function [smEnTC] = background(smEnTC)
            
            bgfn = 'July12_40Vg_DG80_background.txt';
            bgfID = fopen(bgfn, 'r');
            bgdata = fscanf(bgfID, '%f');
            fclose(bgfID);

            chs = 0:1:4096;

            bgHist = [];

            bgi = 1;
            bgk = 1;
            while bgi <= length(bgdata)
                if nonzeros(bgdata(bgi))
                    bgj = 1;
                    while bgj <= bgdata(bgi)
                        bgHist(bgk) = chs(bgi);
                        bgk = bgk + 1;
                        bgj = bgj + 1;
                    end
                end
                bgi = bgi + 1;
            end

            bgHist = (662 * bgHist) / 675;
%             bgHist 

            bgitr = 1;
            while bgitr <= length(bgHist) / 3.75
                smEnTC(end+1) = bgHist(bgitr);
                bgitr = bgitr + 1;
            end

        end

        
        % Filtering low- and high-energy events
        filterCheck = 1;
        while filterCheck <= length(smEnTC)
            if smEnTC(filterCheck) < 50
                smEnTC(filterCheck) = [];
            else
                filterCheck = filterCheck + 1;
            end
        end
        
        filterCheck = 1;
        while filterCheck <= length(smEnTC)
            if smEnTC(filterCheck) > 800
                smEnTC(filterCheck) = [];
            else
                filterCheck = filterCheck + 1;
            end
        end
        
        histogram(smEnTC, 800)
        title('Energy spectrum with smeared peak, binned');
        xlabel('Energy, keV');
        ylabel('Counts');
        grid on;
    end

end