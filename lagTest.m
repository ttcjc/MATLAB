%%

clear variables;
close all;
clc;

%%

caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/run/Run_Test';
% caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_2022/Windsor_SB_wW_Upstream_SC';

namePos = max(strfind(caseFolder, '/')) + 1;
caseName = caseFolder(namePos:end);

cloudName = 'kinematicCloud';

plane = false;
surface = false;
volume = true;

nProc = 4;

%%

[timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');

clc;

%%

[LagProps, LagDataPlane, LagDataSurface, LagDataVolume] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                                            plane, surface, volume, ...
                                                                            timeDirs, deltaT, timePrecision, nProc);

%%

clearvars -except particleProps particleDataVolumetric