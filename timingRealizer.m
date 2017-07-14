function [coinTr, TCen] = timingRealizer()
    
    % NOTE: File read-in currently does not work with Octave

    % Constants for decay (Cs-137)
    % Assuming source strength of 0.11 microCi
    lam = 0.693 / (30 * 3.154 * (10^7));
    N = (0.11 * 3.7*(10^4)) / lam;
    act = lam * N;

    % Grabbing data from G4 PRISM_Sim output file, transforming into matrix
    % pr = 'Enter file name.';
    % fn = input(pr, 's');
    fn = 'output_662keV_1Det_cone_phi0_theta90_abgd.txt';
    fID = fopen(fn, 'r');
    line = fgetl(fID)
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
    TCen = en;
    
    % Coincidence considerations
    i = 1;
    ii = 1;
    
    coinTr = evNum;
    timeAdj = time;
    timeSaver = 0;
    
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
    
    % Timing adjustment for lack of true "global" time in simulation
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

    noTrTCen = nonzeros(TCen);
    
%     histogram(noTrTCen, 1024)


%     helpVec = TCen;
%     itn = 1;
%     while itn <= length(helpVec)
%         if helpVec(itn) ~= 0
%             helpVec(itn) = 1;
%         end
%         itn = itn + 1;
%     end
%     noTrevNum = find(helpVec);
%     noTrhitNum = find(helpVec);
%     noTrtrackID = find(helpVec);
%     noTrdetID = find(helpVec);
%     noTrprocess = find(helpVec);
%     noTrDOI = find(helpVec);
%     noTrPhi = find(helpVec);
%     noTrTheta = find(helpVec);
%     noTrHP = find(helpVec);
%     noTrtime = find(helpVec);
%     
%     noTrdata = [noTrevNum, noTrhitNum, noTrtrackID, noTrTCen, noTrdetID, noTrprocess, noTrDOI, noTrPhi, noTrTheta, noTrHP, noTrtime];
%     
%     newfn = [fn(1:end-4) '_revTC' fn(end-3:end)];
%     newfID = fopen(newfn, 'wt');
%     fprintf(newfID, '%s', line);
%     for ifile = 1:size(noTrdata, 1)
%         fprintf(newfID, '%f', noTrdata(ifile, 1:end));
%         fprintf(newfID, '\n');
%     end
%     fclose(newfID);
    
end