%% Preamble
% ----
% General Purpose Script Preamble

evalc('delete(gcp(''nocreate''));');

clear variables;
close all;
clc;

if exist('/mnt/Processing/Data', 'dir')
    saveLoc = '/mnt/Processing/Data';
else
    saveLoc = '~/Data';
end

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures