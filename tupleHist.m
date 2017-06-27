function [coinAll, coinAllEn, coinTr, coinEn, en, noNoiseEn, evNum, timeAdj] = tupleHist()
    
    % Grabbing data from G4 PRISM_Sim output file, transforming into matrix
    % pr = 'Enter file name.';
    % fn = input(pr, 's');
    fn = 'output_662keV_rand_cone_SA1D_HP912.txt';
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
    HP = dataNC(2:end, 8);
    time = dataNC(2:end, 9);
    
    % Adjustments for "noise" from sim - largely result of ...
    % unknown scatterings in G4
    % Histograms - more can be implemented for more of the data available
    noNoiseEn = nonzeros(en);
%     histogram(noNoiseEn, 256)
%     title('Energy, Binned');
%     xlabel('Energy, keV');
%     ylabel('Counts');
%     grid on;
%     fprintf('Press any key to continue.\n');
%     pause
    tabbedEn = tabulate(round(noNoiseEn));
    ROI;
    fprintf('Press any key to continue.\n');
    pause
    histogram(DOI, 10)
    title('Depth of Interaction, Binned');
    xlabel('Depth of Interaction (mm)');
    ylabel('Counts');
    grid on;
    
    % Grabbing the ROI from the "No-Noise" energy spectrum and ...
    % performing some Gaussian smearing as well as a fit
    function ROI()
        [peakCts, peakInd] = max(tabbedEn(:, 2));
        peakEn = tabbedEn(peakInd, 1);
        n = 1;
        while n <= peakCts
            vecROI(n) = peakEn;
            n = n + 1;
        end
        nFilter = [];
        while length(nFilter) ~= length(vecROI)
            filterFiller = randn / 128; 
            % tightening spread of randn takes resolution at FWHM to  ...
            % approximately 2-percent (13.24 keV is exactly 2%)
            
            if filterFiller
                nFilter = [nFilter, filterFiller + 1];
            end
        end
        smVecROI = nFilter .* vecROI;
        histfit(smVecROI, round(length(unique(smVecROI)) / 16))
        title('Fitted peak with randomly-selected normal smearing');
        xlabel('Energy, keV');
        ylabel('Counts');
        grid on;
        fprintf('Press any key to continue.\n');
        pause
        histSmeared(peakEn, smVecROI);
        resp = input('Is there another peak? (y / n)\n', 's');
        if strcmp(resp, 'y')
            try
                j1 = 1;
                while j1 <= 9
                    tabbedEn(peakInd + j1 - 5, 2) = 0;
                    j1 = j1 + 1;
                end
            end
            ROI;
        elseif strcmp(resp, 'n')
            fprintf('Moving on ...\n');
        end
    end

    % Integrating smeared peak back into full spectrum
    function histSmeared(peakEn, smVecROI)
        x = 1;
        y = 1;
        smEn = noNoiseEn;
        while x <= length(noNoiseEn)
            if smEn(x) == peakEn
                smEn(x) = smVecROI(y);
                y = y + 1;
            end
            x = x + 1;
        end
        histogram(smEn, 256)
        title('Energy spectrum with smeared peak, binned');
        xlabel('Energy, keV');
        ylabel('Counts');
        grid on;
        fprintf('Press any key to continue.\n');
        pause
    end
    
    % Coincidence considerations
    i = 1;
    coinTr = evNum;
    timeAdj = time;
    while i < length(evNum)
        if (evNum(i) == evNum(i+1)) && (detID(i) == detID(i+1))
            coinTr(i) = 1;
            coinTr(i+1) = 1;
            timePl = timeAdjUq(i);
            timeAdjNonUq(i+1, timePl);
            i = i + 2;
        elseif (evNum(i) == evNum(i-1)) && (detID(i) == detID(i-1))
            coinTr(i) = 1;
            timeAdjNonUq(i, timePl);
            i = i + 1;
        else
            coinTr(i) = 0;
            timeAdjUq(i);
            i = i + 1;
        end
    end
    if coinTr(end) ~= 1
        coinTr(end) = 0;
        timeAdjUq(i)
    end
    
    coinEn = coinTr .* en;
    coinEn = coinEn(find(coinTr));
    
    % Timing adjustment for lack of true "global" time in simulation
    function timeAdjNonUq(ind, timePl)
        timeAdj(ind) = timeAdj(ind) + timePl;
    end
    
    % consider that 3220 gammas are emitted by one gram of Cs-137 in 1 ns
    function [timePl] = timeAdjUq(ind)
        if evNum(ind) < 3220
            timePl = rand;
            timeAdj(ind) = timeAdj(ind) + timePl;
        else
            timeAdd = floor(evNum(ind) / 3220);
            timePl = rand + timeAdd;
            timeAdj(ind) = timeAdj(ind) + timePl;
        end
    end
    
    % identifying "false" coincidences based on time (including ...
    % adjustments)
    k = 1;
    coinAll = evNum;
    while k < length(timeAdj)
        if ((timeAdj(k) + 0.2) >= timeAdj(k+1)) && ((timeAdj(k) - 0.2) <= timeAdj(k+1)) && (detID(k) == detID(k+1))
            coinAll(k) = 1;
            coinAll(k+1) = 1;
            k = k + 2;
        elseif k > 1
            if ((timeAdj(k) + 0.2) >= timeAdj(k-1)) && ((timeAdj(k) - 0.2) <= timeAdj(k-1)) && (detID(k) == detID(k-1))
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
    
%     histogram(timeAdj, 100)
%     title('Time, Binned');
%     xlabel('Time, ns');
%     ylabel('Counts');
%     grid on;
end