function tupleHist()
    
    % NOTE: File read-in currently does not work with Octave

%     profile on;

    % Constants for decay (Cs-137)
    % Assuming source strength of 10 microCi
    lam = 0.693 / (30 * 3.154 * (10^7));
    N = (10.08 * 3.7*(10^4)) / lam;
    act = lam * N / 1000000000;

    % Grabbing data from G4 PRISM_Sim output file, transforming into matrix
    % pr = 'Enter file name.';
    % fn = input(pr, 's');
    
    % Code block designed for reading Geant4 output files and may ...
    %   require slight tweaking for different formats
    fn = 'output_662keV_1Det_isotropic_phi360_theta90.txt';
    fID = fopen(fn, 'r');
    line1 = fgetl(fID);
    fclose(fID);
    isStrCol = isnan(str2double(regexp(line1, '[^\t]+', 'match')));
    
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
    
    % Uncomment below block for a non-processed spectrum 
%     nonZen = nonzeros(en);
%     
%     [Nde, Xde] = hist(nonZen, 1024);
%     Hde = bar(Xde, Nde, 1);
%     Hde.FaceAlpha = 0.0;
% 
%     Nde = conv(Nde, ones(1, 1), 'same') / 1;
%     Hde = line(Xde, Nde);
%     set(Hde, 'color', 'g', 'linewidth', 0.8)
% 
%     title('Energy spectrum, binned (Simulation Data)');
%     xlabel('Energy, keV');
%     xlim([0 1024]);
%     ylabel('Counts');
%     grid on;
%     fprintf('Press any key to continue.');
%     pause;
    
    % Uncomment for DOI histogram
%     histogram(DOI, 10)
%     title('Depth of Interaction, Binned');
%     xlabel('Depth of Interaction (mm)');
%     ylabel('Counts');
%     grid on;
%     fprintf('Press any key to continue.');
%     pause;
    
    % Coincidence Considerations ------------------------------------------
    % Initializing indexing variables and modified arrays
    i = 1;
    ii = 1;
    
    coinTr = evNum;
    timeAdj = time;
    timeSaver = 0;
    TCen = en;
    
    % Considering events that miss the detector geometry, adjusting ...
    %   timing accordingly
    while ii <= evNum(end) && i <= length(evNum)
        % Below block is legacy code - replaced with much faster ...
        % indexing farther below
%         if ~any(i==evNum)
%             while ~any(ii==evNum)
%                  timeSaver = timeAdjGhost(timeSaver);
%                  ii = ii + 1;
%             end
%         end
        
        if evNum(i) == evNum(1)
            evInc = evNum(i+1) - 1;
                ghostEvs = evNum(i+1) - 1;
                while ghostEvs > 0
                    timeSaver = timeAdjGhost(timeSaver);
                    ghostEvs = ghostEvs - 1;
                end
            ii = ii + evInc;
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
                
                timeAdj(i) = timeAdj(i-1) + tTracker;
                i = i + 1;
                ctr = ctr + 1;
            end
            
        end
        
        if evNum(i-1) ~= evNum(end)
            if evNum(i-1) ~= evNum(i)
                evInc = evNum(i) - evNum(i-1);
                ghostEvs = evNum(i) - evNum(i-1) - 1;
                while ghostEvs > 0
                    timeSaver = timeAdjGhost(timeSaver);
                    ghostEvs = ghostEvs - 1;
                end
                ii = ii + evInc;
            end
        end
        
        if i ~= 1
            ctr = 1;
            try
                while evNum(i-1) == evNum(i)
                    tTracker = time(i) - time(i-1);

                    coinTr(i) = evNum(i);
                    coinTr(i-1) = evNum(i-1);

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
                        
                        timeAdj(i) = timeAdj(i-1) + tTracker;
                        i = i + 1;
                        ctr = ctr + 1;
                    end
                end
            end
        end
        ii = ii + 1;
        
    end
    
    % Primarily for debugging purposes
    coinEn = coinTr .* en;
    coinEn = coinEn(find(coinTr));
    
    % Timing adjustment functions for lack of true "global" time in ...
    % Geant4 simulation readout
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

    % Identifying "false" coincidences based on time (including ...
    %   adjustments)
    % Also adjusting energy vector due to limitations of real detector ...
    %   (summing multi-track events contained within detector)
    
    k = 1;
    coinAll = evNum;
    coinAllID = 1;
    while k < length(timeAdj)
        % Note: 1.3 mircoseconds discussed in meeting as being current
        % Goal = 0.45 microSec window
        
        % Units of time are ns, adjust number to vary window for chance ...
        %   coincidence
        if ((timeAdj(k) + 1000) >= timeAdj(k+1)) && ((timeAdj(k) - 1000) <= timeAdj(k+1))
            
            coinAll(k) = coinAllID;
            coinAll(k+1) = coinAllID;
            
            % Finding a "weighted" DOI for multi-deposition events
            holdEn1 = TCen(k);
            holdEn2 = TCen(k+1);            
            enSum = holdEn1 + holdEn2;
            enAr = [holdEn1, holdEn2];
                        
            ctr1 = 2;
            try
                while ((timeAdj(k) + 1000) >= timeAdj(k+ctr1)) && ((timeAdj(k) - 1000) <= timeAdj(k+ctr1))
                    coinAll(k+ctr1) = coinAllID;
                    fracEnAr = [];
                    
                    holdEnX = TCen(k+ctr1);
                    enSum = enSum + holdEnX;
                    
                    enAr(end+1) = holdEnX;
                    
                    ctr1 = ctr1 + 1;
                end
            end
            if enSum ~= 0
                fracEnAr = enAr / enSum;
            else
                fracEnAr = enAr;
            end

            doiWeighted = 0;
            fracCtr = 1;
            while fracCtr <= length(fracEnAr)
                doiHold = DOI(k+fracCtr-1);
                doiFracHold = fracEnAr(fracCtr) * doiHold;
                doiWeighted = doiWeighted + doiFracHold;

                fracCtr = fracCtr + 1;
            end
            
            % Debugging blocks            
            if doiWeighted > 10
                fprintf('something gone wrong (more than 10)\n');
                saveAr = fracEnAr;
            elseif (doiWeighted < 0) && (doiWeighted ~= -1)
                fprintf('something gone wrong (less than 0)\n');
                saveAr = fracEnAr;
            end
            if isnan(doiWeighted)
                fprintf('something gone wrong (NaN)\n');
                saveAr = fracEnAr;
            end
            
            % Modifying energy and DOI vectors
            TCen(k) = TCen(k) + TCen(k+1);
            TCen(k+1) = 0;
            
            ctr11 = ctr1;
            ctr12 = ctr1;
            try
                while ctr11 > 2
                    TCen(k) = TCen(k) + TCen(k+ctr11-1);
                    TCen(k+ctr11-1) = 0;
                    ctr11 = ctr11 - 1;
                end
            end
            
            DOI(k+1) = -1;
            
            try
                while ctr12 > 2
                    DOI(k+ctr12-1) = -1;
                    ctr12 = ctr12 - 1;
                end
            end
            
            DOI(k) = doiWeighted;
            
            k = k + ctr1;
            coinAllID = coinAllID + 1;
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
    
    % Un-comment block for debugging timing adjustment
%     histogram(timeAdj, 100)
%     title('Time, Binned');
%     xlabel('Time, ns');
%     ylabel('Counts');
%     grid on;
%     fprintf('Press any key to continue.\n');
%     pause
    
    % Adjusting data vectors ----------------------------------------------
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
    
    % If everything went correctly with energy summing, noTrDOI should ...
    %   have NO values less than 0 - uncomment below statement for ...
    %   debugging
    noTrDOI = DOI(find(helpVec));
    
%     if any(noTrDOI < 0)
%         fprintf('Bad news, negatives in en-summing modified DOI vec\n');
%     end
    
    noTrPhi = Phi(find(helpVec));
    noTrTheta = Theta(find(helpVec));
    noTrHP = HP(find(helpVec));
    noTrtime = time(find(helpVec));
    
    % Tabulating energy vector to pick out peak
    tabbedEn = tabulate(round(noTrTCen));
    
    % Processing peak through smearing algorithms -------------------------
    [noiseFWHM, howTight] = ROITC;
    
    function [noiseFWHM, howTight] = ROITC()
        
        % Grabbing peak statistics
        [peakCts, peakInd] = max(tabbedEn(:, 2));
        peakEn = tabbedEn(peakInd, 1);
        n = 1;
        
        % Assuming a measurement of 10 keV-FWHM for electronic noise
        elecNoiseFWHM = 10;
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
        %   target (the noise fwhm)
        function [nFilter, smVecROI, howTight] = tightenerTC(nFilter, howTight)
            nFilter = [];
            while length(nFilter) ~= length(vecROI)
                
                filterFiller = randn / howTight; 
                % Tightening spread of randn takes resolution at FWHM ...
                %   to approximately 2-percent (13.24 keV is exactly ...
                %   2% for Cs-137) initially, will be revised in ...
                %   recursion later

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
            
            % Two while loops to grab energies of the left and right ...
            %   edges of the fwhm
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
            %   than that dictated by the noise
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
        
        % Demonstrating normally-smeared peak region
%         histfit(smVecROI, round(length(unique(smVecROI)) / 16))
%         title('Fitted peak with randomly-selected normal smearing');
%         xlabel('Energy, keV');
%         ylabel('Counts');
%         grid on;
%         fprintf('Press any key to continue.\n');
%         pause
        
        histSmearedTC(peakEn, smVecROI, howTight, elecNoiseSTDDEV);
        
        % Uncomment below code if more than one peak - however, with ...
        %   G4 gun limitations, unsure if this functionality would ever ...
        %   be usable
        % NOTE: Has not been extensively tested
        
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
        
        % Integrating "smeared" energies back into full spect in place ...
        %   monoenergetic peak energies
        % Also smearing of non-peak region
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
            
            % NOTE: Algorithms constructed with consideration of CPG ...
            %   and near-perfect DG optimization
            % NOTE: Assuming detector positioning of cathode above ...
            %   anode (if center of crystal is origin, cathode is at ...
            %   +5mm on z-axis, anode is at -5mm on z-axis)
            
            % Numbers based on Michelle's dissertation
            muTaoEl = 0.01;
            muTaoHo = 0.0008;
            
            bias = 1000;
            
            detdepth = 1;
            
            % Modify this variable for different isotope backscatter ...
            %   peak energies
            backScatPeakEn = 180;
            
            clct = 1;
            dbgEnVec = [];
            dbgOrigEnVec = [];
            dbgDOIVec = [];
            dbgEnVec2 = [];
            
            while clct <= length(smEnTC)
                
                % Using Hecht's equation:
                CCEmax = ((muTaoEl * bias)/(detdepth) * (1 - exp(-(detdepth - (0 * 0.1)) / (muTaoEl * bias)))) + ((muTaoHo * bias)/(detdepth) * (1 - exp(-(0 * 0.1) / (muTaoHo * bias))));
                
                enCL = smEnTC(clct) * CCEmax;

                % Adjusting for imperfect DG setting - linearly ...
                %   electron trapping with depth, up to 0.5 percent at ...
                %   DOI = 10mm
                
                enMod1 = 1 - (noTrDOI(clct) * 0.0005);
                enCL = enCL * enMod1;
                
                % Exponential degradation of signal for interactions ...
                %   near anode
                if noTrDOI(clct) <= 1
                    expDeg = noTrDOI(clct) - 1;
                    enMod2 = exp(expDeg);
                    enCL = enCL * enMod2;
                end
                
                if noTrDOI(clct) <= 2
                    if enCL > (peakEn * 0.9516) - 50
                        enCL = enCL * 0.98;
                        dister = randn / 10;
                        while dister > -0.2040
                            dister = randn / 10;
                        end
                        dister = dister + 1;
                        enCL = 135 + (enCL * dister);
                    end
                    
                end
                
                if noTrDOI(clct) > 2 && noTrDOI(clct) < 4
                    if enCL > (peakEn * 0.9516) - 20
                        if rand > 0.90
                            disterLowEn = randn / 6;
                            enCL = (enCL * (60 / smEnTC(clct))) * (disterLowEn + 1);
                        end
                    end
                    
                    if enCL > (backScatPeakEn * 0.9516) - 80 && smEnTC(clct) < (backScatPeakEn * 0.9516) + 20
                        if rand > 0.60
                            disterLowEn = randn / 6;
                            enCL = (enCL * (60 / smEnTC(clct))) * (disterLowEn + 1);
                        end
                    end
                    
                    if enCL > (peakEn * 0.9516) - 5
                        if rand < 0.2
                            enCL = enCL * 0.7251;
                            disterHiTail = randn / 12;
                            while disterHiTail < 0.0860
                                disterHiTail = randn / 12;
                            end
                            disterHiTail = disterHiTail + 1;
                            enCL = (enCL * disterHiTail) - 44;
                        elseif rand > 0.9925
                            enCL = enCL / 100;
                        end
                    end
                    
                end
                
                % Adjusting for possibility for collection of holes in ...
                %   region near to cathode
                if noTrDOI(clct) >= 8
                    
                    if enCL > (peakEn * 0.9516) - 50
                        enCL = enCL * 0.98;
                        dister = randn / 10;
                        while dister > -0.2040
                            dister = randn / 10;
                        end
                        dister = dister + 1;
                        enCL = 135 + (enCL * dister);
                    end
                        
                end
                
                if noTrDOI(clct) > 6 && noTrDOI(clct) < 8
                    if smEnTC(clct) > (peakEn * 0.9516) - 20
                        if rand > 0.95
                            disterLowEn = randn / 6;
                            enCL = (enCL * (60 / smEnTC(clct))) * (disterLowEn + 1);
                        end
                    end
                    
                    if smEnTC(clct) > backScatPeakEn - 80 && smEnTC(clct) < backScatPeakEn + 20
                        if rand > 0.75
                            disterLowEn = randn / 6;
                            enCL = (enCL * (60 / smEnTC(clct))) * (disterLowEn + 1);
                        end
                    end
                    
                    if enCL > (peakEn * 0.9516) - 5
                        if rand < 0.2
                            enCL = enCL * 0.7251;
                            disterHiTail = randn / 12;
                            while disterHiTail < 0.0860
                                disterHiTail = randn / 12;
                            end
                            disterHiTail = disterHiTail + 1;
                            enCL = (enCL * disterHiTail) - 44;
                        elseif rand > 0.985
                            enCL = enCL / 100;
                        end
                    end

                end
                
                % Filtering out excess counts in compton continuum and ...
                %   low-energy peak
                if enCL < 65
                    if rand < 0.25
                        enCL = enCL / 100;
                    end
                end
                if (enCL < 475.8) && (enCL > (240 * 0.9516))
                    enLivMod = (enCL - (220 * 0.9516)) / 1000;
                    if rand < enLivMod
                        enCL = enCL / 100;
                    end
                end
                
                if enCL >= 600 && enCL <= 645 && smEnTC(clct) > 650
                    dbgEnVec(end+1) = enCL;
                    dbgOrigEnVec(end+1) = smEnTC(clct);
                    dbgDOIVec(end+1) = noTrDOI(clct);
                end
                
                smEnTC(clct) = enCL;
                clct = clct + 1;
            end
            
        end
        
        % Uniformly adjusting for charge loss shift
        smEnTC = smEnTC * 1.0509;
        
        % Adding the measured background into the spectrum
        function [smEnTC] = background(smEnTC)
            
            bgfn = 'July12_40Vg_DG80_background_tenmin.txt';
            bgfID = fopen(bgfn, 'r');
            bgdata = fscanf(bgfID, '%f');
            fclose(bgfID);
            
            % Change this line if background taken with different ...
            %   number of channels
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
            
            % Change this line for different isotopes/calibrations
            bgHist = (662 * bgHist) / 675;

            bgitr = 1;
            while bgitr <= length(bgHist) / 3.75
                smEnTC(end+1) = bgHist(bgitr);
                bgitr = bgitr + 1;
            end

        end
        
        % Filtering low- and high-energy events for better comparison ...
        %   to filtered experimental data
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
            if smEnTC(filterCheck) > 1024
                smEnTC(filterCheck) = [];
            else
                filterCheck = filterCheck + 1;
            end
        end
        
        % Uncomment this line and comment lines 664 thru 670 for a bar hist
        histogram(smEnTC, 1024)
        
%         [Nd, Xd] = hist(smEnTC, 1024);
%         Hd = bar(Xd, Nd, 1);
%         Hd.FaceAlpha = 0.0;
% 
%         Nd = conv(Nd, ones(1, 1), 'same') / 1;
%         Hd = line(Xd, Nd);
%         set(Hd, 'color', 'g', 'linewidth', 1.25)
        
        title('Energy spectrum with smeared peak, binned (Simulation Data)');
        xlabel('Energy, keV');
        xlim([0 1024]);
        ylabel('Counts');
        grid on;
        
        % Uncomment these two lines if you wish to use the profile viewer
%         fprintf('Press any key to continue.\n');
%         pause
    end

%     profile viewer

end