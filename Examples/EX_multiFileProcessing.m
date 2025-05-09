% Demonstrating how to load many .mat files and perform arbitrary processing on them 
% For this example, get all peaks for all mean sweep regions in all files
% Modify the code inside the loop of 'i' to process all data
clc
clear
close all
%% COUNT FILES AND INTIALIZE VARIABLES
% How to split a string variable over multiple lines:
%{
dataFolder = ['D:\Google Drive\Career\Career - MAIN\2015-2019\' ... 
'Dr bickford laboratory\Work\Projects\Clampex data analysis\code\Examples\ExData\'];
%}
dataFolder = uigetdir();                       % Select the folder 'ExData' for this example
if dataFolder == 0; clear; return; end         % If user cancels folder prompt, close without throwing an error
numFiles = length(dir([dataFolder '\*.mat'])); % Find number of .mat files in chosen folder
dir = dir(dataFolder);                         % Used in loop to get file names for loading
myPeaks(1,:) = 0;                              % Initialize a variable to store your data in
%% LOOP OVER FILES
for fileID = 1:1:numFiles % Loop over files
    currFile = dir(fileID+2,1).name;         % Gets the first file (ignores two system files '.' and '..')
    load(currFile);
    %% YOUR PER FILE PROCESSING HERE
    allMeanStructs = sweepData(:, end);      % Uses end since the mean is the last column. (row, column)
    for i = 1:1:numel(allMeanStructs)        % Loop over structs
        currRegStr = allMeanStructs(i);      % Current region struct
        myPeaks(i, fileID) = currRegStr.pkV; % Each column of myPeaks will be the peaks from the mean sweep a single file        
    end
end
%% SAVE TO FILE (IF DESIRED)
matName = 'exampleData';
save(matName, 'myPeaks'); % Save myPeaks variable to a file called exampleData in the current folder
disp('DONE');