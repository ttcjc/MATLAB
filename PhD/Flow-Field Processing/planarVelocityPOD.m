%% Planar Velocity POD Calculator v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===================================');
disp ('Planar Pressure POD Calculator v1.0');
disp ('===================================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialisation

[caseName, probeData, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(true, normalise, nProc);    

disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)