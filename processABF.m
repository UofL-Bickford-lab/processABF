%{
% processABF.m 
% PURPOSE: Take data from a series of .abf files (clampex) and analyze the
% data relative to laser regions.
%
% INPUTS: Folder containing .abf files (from user prompt). Also see
% 'CONTROLS' section in code.
%
% OUTPUTS: The final data in the form of structs in a file (.mat).
%
% DEPENDENCIES:  Basic MATLAB install (R2019b or later due to tiledlayout function).
% Signal processing toolbox for peakfinding function.
% Also function abf2load.m (abfload.m is old version, won't work with newer
% versions of clampex abf files).
%
% AUTHOR: David C Alston (dalston2428@gmail.com) 2020
%
% NOTES:
%   - Assumes ~5/0V laser pulses (looks for diff(laserVoltage)>1 to find edges)
%   - Puts the output files in the current folder for Matlab.
%   - Each column of sweepData represents one sweep. The last column is analysis of the average sweep.
%       -- So if you select sweeps 1 and 3 from three sweeps, the sweepData
%          will have three columns (data for sweep 1, data for sweep 3, and data for the mean sweep)
%   - Each single element of sweepData represents one laser region + the specified search window and baseline.
%       -- .window is the data searchWidthms seconds after the rising edge of the laser for this laser region.
%       -- .baseline is the voltage/time in seconds before the laser region (see CONTROLS).
%       -- .baselineAvg is the mean of the .baseline voltage.
%       -- .lasStartTimeABS is the time of the laser rising edge relative to the start of the data
%       -- .pkVTimeABS is the time of the peak in seconds relative to the start of the data
%       -- .pkVTimeREL is the time of the peak in seconds relative to the laser rising edge
%   - For origData.data, it is of the format (data, channelNum, sweepNum)
%       -- For example, myVar = origData.data(:,:,1) is all the data from sweep one.
%       -- myVar = origData.data(:,2,1) is all the data from channel two in sweep one.
%}
clc
clear
close all
%% CONTROLS
searchWidthms = 30;                % milliseconds after start of a laser pulse
baselineWidthms = 2;               % milliseconds before start of laser pulse (not including pulse start index)
lasChanName = 'IN 4';         % Channel containing laser data
spikeChanName = 'IN 0';       % Channel containing voltage data
minOrMax = 'min';                  % 'min' = find troughs, 'max' = find peaks
%% MAIN PROGRAM
basePath = uigetdir([], 'Select folder containing .abf files');
%basePath = '\Your\Path\To\Files\'; % For debugging
if basePath == 0; clear; return; end % If user cancels folder prompt, close without throwing an error
numFiles = length(dir([basePath '\*.abf']));
if numFiles == 0; beep; disp('ERROR:: No .abf files found in folder. Closing...'); clear; return; end
dir = dir(basePath);
for i = 1:1:numFiles
    %% Check inputs from CONTROLS and load raw data from file
    clc
    clearvars -except i basePath numFiles dir searchWidthms baselineWidthms lasChanName spikeChanName minOrMax
    currFile = dir(i+2,1).name;
    [origData.data, origData.sampIntmicroSec, origData.fileInfo] = abf2load(strcat(basePath, '\', currFile));
    searchWidthInd = round((searchWidthms/1E3) / (origData.sampIntmicroSec/1E6));
    baselineWidthInd = round((baselineWidthms/1E3) / (origData.sampIntmicroSec/1E6));
    if searchWidthInd < 2
        clc
        clear
        beep;
        disp('ERROR:: Search window width would be less than two data points. Check your chosen search width. Closing...');
        return
    end
    if baselineWidthInd < 2
        clc
        clear
        beep;
        disp('ERROR:: Baseline width would be less than two data points. Check your chosen baseline width. Closing...');
        return
    end    
    numSweeps = size(origData.data, 3);
    numChannels = size(origData.data, 2);
    %% Find Channels
    spikeFlag = 0;
    lasFlag = 0;
    for chanNum = 1:1:numChannels % Find laser and spike channels and grab data
        chanName = char(origData.fileInfo.recChNames(chanNum));
        if contains(chanName, lasChanName, 'IgnoreCase', true)
            lasChanIndx = chanNum;
            lasFlag = 1;
        end
        if contains(chanName, spikeChanName, 'IgnoreCase', true)
            spikeChanIndx = chanNum;
            spikeFlag = 1;
        end
        if (spikeFlag && lasFlag); break; end % Found both channels so done (break loop)
    end
    if (~lasFlag || ~spikeFlag)
        disp('ERROR:: Laser or spike channel name not found for current file. Check your channel names under program controls. Channel names found:');
        for chanNum = 1:1:numChannels
            chanName = char(origData.fileInfo.recChNames(chanNum));
            disp(chanName);
        end
        beep;
        clear
        return
    end
    %% Display and select sweeps for averaging/analysis
    allSweeps = zeros(size(origData.data, 1), numSweeps);
    fig = tiledlayout(1, numSweeps);
    for n = 1:1:numSweeps
        allSweeps(:,n) = origData.data(:, spikeChanIndx, n);
        currPlot = nexttile;
        plot(allSweeps(:,n));
        title(currPlot, int2str(n));
        xlabel('Index');
        ylabel('mV');
    end
    plotFig = gcf;
    plotFig.WindowState = 'maximized';
    clearvars list
    list{:} = '0';
    for b = 1:1:numSweeps; list{b} = strcat('Sweep_', int2str(b)); end
    [indx, validEntry] = listdlg('ListString', list);
    if ~validEntry; disp('INFO:: User cancelled sweep selection. Closing...'); close all; clear; return; end
    close
    selSweeps = zeros(size(origData.data, 1), numel(indx)); % Selected sweeps
    for m = 1:1:numel(indx)
        currIndx = indx(m);
        selSweeps(:, m) = allSweeps(:, currIndx);
    end
    selSweeps(:,m+1) = mean(selSweeps, 2);
    meanSweep(:,1) = selSweeps(:,m+1); % For saving to .mat later on
    %% Analyze chosen sweeps + the mean sweep
    laserData = origData.data(:,lasChanIndx, 1); % Assuming laser is same across sweeps
    for sweepNum = 1:1:numel(indx)+1 % Analyze chosen sweeps + mean sweep (mean of chosen sweeps)
        currVoltage = selSweeps(:, sweepNum);
        clearvars timeArrSec
        timeArrSec(:,1) = zeros(numel(laserData), 1);
        for n = 2:1:numel(laserData); timeArrSec(n,1) = timeArrSec(n-1,1)+(origData.sampIntmicroSec/1E6); end
        dy = abs(gradient(laserData));
        [~, lasIndRange] = findpeaks(dy, 'MinPeakDistance',2, 'MinPeakHeight',1); % Requires signal processing toolbox
        numRegions = numel(lasIndRange)/2;
        count = 0;
        for l = 2:2:numel(lasIndRange)
            count = count+1;
            lasInd(1, count) = lasIndRange(l-1);
            lasInd(2, count) = lasIndRange(l);            
        end
        lasStruct(1).regionInd = 0;
        for c = 1:1:size(lasInd, 2)
            startInd = lasInd(1, c);
            endInd = lasInd(2, c);
            numPts = endInd - startInd;
            for v = 1:1:numPts+1
                lasStruct(c).regionInd(v, 1) = startInd+v-1;
            end
        end
        baselineAvgs = zeros(numRegions, 1);
        for m = 1:1:numRegions % Get baseline just before each laser region
            lasStartInd = lasStruct(m).regionInd(1);
            baselineStartInd = lasStartInd - baselineWidthInd - 1; % Want to end 1 index before index where laser starts
            for b = 1:1:baselineWidthInd
                baselineInds(b, m) = baselineStartInd + b; % Each column is the selected time before the laser (for baseline)
            end       
            baselineAvgs(m, 1) = mean(currVoltage(baselineInds(:, m)));
        end
        baselineVoltages = currVoltage(baselineInds);
        for b = 1:1:numRegions % Get and process the window after each laser region
            startInd = lasStruct(b).regionInd(1); % Start search at rising edge of laser
            endInd = startInd+searchWidthInd-1;
            lasRisingEdge = timeArrSec(startInd);
            inds = startInd:endInd;
            voltage = currVoltage(inds);
            if strcmp(minOrMax, 'max')
                [pkVal, K] = max(voltage);
            else
                [pkVal, K] = min(voltage);
            end
            sweepData(b, sweepNum).peakMode = minOrMax;
            sweepData(b, sweepNum).window(:,1) = voltage; %#ok<*SAGROW>
            sweepData(b, sweepNum).window(:,2) = timeArrSec(inds);
            sweepData(b, sweepNum).pkV = pkVal;
            sweepData(b, sweepNum).pkVTimeABS = timeArrSec(inds(K));
            sweepData(b, sweepNum).pkVTimeREL = sweepData(b, sweepNum).pkVTimeABS - lasRisingEdge;
            sweepData(b, sweepNum).lasStartTimeABS = lasRisingEdge;
            sweepData(b, sweepNum).baseline(:, 1) = baselineVoltages(:, b);
            sweepData(b, sweepNum).baseline(:, 2) = timeArrSec(baselineInds(:, b));
            sweepData(b, sweepNum).baselineAvg = baselineAvgs(b);
        end
    end
    meanSweep(:,2) = timeArrSec;
    %% Save to file
    chosenSweeps = indx;
    outName = erase(currFile, '.abf');
    save(outName, 'meanSweep', 'chosenSweeps', 'origData', 'sweepData', 'searchWidthms', 'baselineWidthms', 'lasChanName', 'spikeChanName');
end
clear
disp('Finished analyzing .abf files');