% Demonstrating how to load a single file and get the half width in a variable
clc
clear
close all
%% CONTROLS
addpath('ExData'); % Since the ExData folder is in the current folder this is all that is needed
myFile = '2019_10_18_0003.mat';
halfWidthY = -66.6; % User selected Y value (mV)
%% LOAD FILE
load(myFile); % Load the .mat file into the workspace
clearvars dataFolder myFile % Just to declutter the workspace a little bit
%% MAIN PROGRAM
% Using first laser pulse and mean sweep:
region1V = sweepData(1, 3).window(:,1); % Just grab voltage data
[X0] = intersect1D_F(region1V, halfWidthY); % All X intersect indices with chosen Y value
if isempty(X0)
    beep;
    disp('No intersections found. Try changing your halfWidthY variable. Closing...');
    return
end
firstIntersect = X0(1);
for n = 2:1:numel(X0) % Find the next intersect that is > 2 indicies away from firstIntersect
    currDiff = X0(n) - firstIntersect;
    if currDiff > 2 % Found it, break out of the for loop
        finalDelta = currDiff; % Final half width
        break;
    end
end
%% PRINTING TO WINDOW AND PLOTTING
fprintf('The final half width found was %i indices wide \n', finalDelta);
plot(region1V);
hold on
yline(halfWidthY, '-.r');
hold on
for i = 1:1:numel(X0)
    hold on
    xline(X0(i));   
end