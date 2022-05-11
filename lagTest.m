clear variables;
close all;
clc;

%%

% caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/run/Run_Test';
caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_2022/Windsor_SB_wW_Upstream_SC';

%%

cloudName = 'kinematicCloud';

%%

[timeDirs, ~] = timeDirectories(caseFolder, 'global');

%%

% [particleData, ~] = readLagrangianDataPlanar(caseFolder, timeDirs, cloudName, 32);
% [particleData, ~] = readLagrangianDataPlanar(caseFolder, timeDirs, cloudName, 400);

%%

[particleData, ~] = readLagrangianDataVolume(caseFolder, timeDirs, cloudName, 4);