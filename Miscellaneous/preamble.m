%% Preamble
% ----
% General Purpose Script Preamble

evalc('delete(gcp(''nocreate''));');

clear variables;
close all;
clc;

% Start of Save Location File Path
if exist('/mnt/Processing/Data', 'dir')
    saveLoc = '/mnt/Processing/Data';
else
    saveLoc = '~/Data';
end

geoLoc = '~/CAD/CFD Geometries'; % Start of Geometry File Path

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures