function [smEnTC] = tupleHist()
    
    % NOTE: File read-in currently does not work with Octave

    % Constants for decay (Cs-137)
    % Assuming source strength of 0.10 microCi
    lam = 0.693 / (30 * 3.154 * (10^7));
    N = (0.1008 * 3.7*(10^4)) / lam;
    act = lam * N;

    % Grabbing data from G4 PRISM_Sim output file, transforming into matrix
    % pr = 'Enter file name.';
    % fn = input(pr, 's');
    fn = 'output_662keV_1Det_cone_phi0_theta90.txt';
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
                
                TCen(i-ctr) = TCen(i-ctr) + TCen(i);
                TCen(i) = 0;
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

                    TCen(i-ctr) = TCen(i-ctr) + TCen(i);
                    TCen(i) = 0;

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
                        
                        TCen(i-ctr) = TCen(i-ctr) + TCen(i);
                        TCen(i) = 0;
                        
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
    k = 1;
    coinAll = evNum;
    while k < length(timeAdj)
        if ((timeAdj(k) + 0.005) >= timeAdj(k+1)) && ((timeAdj(k) - 0.005) <= timeAdj(k+1)) && (detID(k) == detID(k+1))
            coinAll(k) = 1;
            coinAll(k+1) = 1;
            k = k + 2;
        elseif k > 1
            if ((timeAdj(k) + 0.005) >= timeAdj(k-1)) && ((timeAdj(k) - 0.005) <= timeAdj(k-1)) && (detID(k) == detID(k-1))
                coinAll(k) = 1;
                k = k + 1;
            else
                coinAll(k) = 0;
                k = k + 1;
            end
        else
            coinAll(k) = 0;
            k = k + 1;
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
    
    helpVec = TCen;
    itn = 1;
    while itn <= length(helpVec)
        if helpVec(itn) ~= 0
            helpVec(itn) = 1;
        end
        itn = itn + 1;
    end
    noTrevNum = find(helpVec);
    noTrhitNum = find(helpVec);
    noTrtrackID = find(helpVec);
    noTrdetID = find(helpVec);
    noTrprocess = find(helpVec);
    noTrDOI = find(helpVec);
    noTrPhi = find(helpVec);
    noTrTheta = find(helpVec);
    noTrHP = find(helpVec);
    noTrtime = find(helpVec);
    
    tabbedEn = tabulate(round(noTrTCen));

    [noiseFWHM, howTight] = ROITC;
    
    function [noiseFWHM, howTight] = ROITC()
        [peakCts, peakInd] = max(tabbedEn(:, 2));
        peakEn = tabbedEn(peakInd, 1);
        n = 1;
        
        % Assuming a measurement of 9.5 keV-FWHM for electronic noise
        elecNoiseFWHM = 9.5;
        elecNoiseSTDDEV = elecNoiseFWHM / 2.355;
        
        % For an average ionization energy in CZT of 5 eV
        qCarrGenN = peakEn / 0.005;
        % For an average pulse height of 0.350 V
        propCons = 0.350 / qCarrGenN;
        statNoiseSTDDEV = propCons * sqrt(qCarrGenN);
        
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
        histSmearedTC(peakEn, smVecROI, howTight);
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
        
    function histSmearedTC(peakEn, smVecROI, howTight)
        x = 1;
        y = 1;
        smEnTC = noTrTCen;
                
        while x <= length(noTrTCen)
            if round(smEnTC(x)) == peakEn
                smEnTC(x) = smVecROI(y);
                y = y + 1;
            else 
                noiseFilter = randn / howTight;
                smEnTC(x) = smEnTC(x) .* (noiseFilter + 1);
            end
            x = x + 1;
        end
        
        smEnTC = background(smEnTC);

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

            bgitr = 1;
            while bgitr <= length(bgHist) / 3.75
                smEnTC(end+1) = bgHist(bgitr);
                bgitr = bgitr + 1;
            end

        end
        
        histogram(smEnTC, 1024)
        title('Energy spectrum with smeared peak, binned');
        xlabel('Energy, keV');
        ylabel('Counts');
        grid on;
    end

end