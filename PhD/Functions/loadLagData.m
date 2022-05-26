%% Independent Lagrangian Data Reader v1.0
% ----
% Executes 'initialiseLagData.m' Without Further Processing
% ----
% Usage: [LagProps, LagDataPlane, LagDataSurface, LagDataVolume] = loadLagData(caseFolder, cloudName, ...
%                                                                              plane, surface, volume);
%        'caseFolder' -> Case Path, Stored as s String
%        'cloudName'  -> OpenFOAM Cloud Name, Stored as a String
%        'plane'      -> Collect Planar Data [True/False]
%        'surface'    -> Collect Surface Data [True/False]
%        'volume'     -> Collect Volume Data [True/False]


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [LagProps, LagDataPlane, LagDataSurface, ...
          LagDataVolume, sampleInterval] = loadLagData(caseFolder, cloudName, ...
                                                       plane, surface, volume)
                                                    
    nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseName = caseFolder(namePos:end);
    
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');
    
    clc;
    
    [LagProps, LagDataPlane, LagDataSurface, ...
     LagDataVolume, sampleInterval] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                        plane, surface, volume, ...
                                                        timeDirs, deltaT, timePrecision, nProc);
end