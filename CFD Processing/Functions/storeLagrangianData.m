%% Lagrangian Data Storage v1.1
% ----
% Executes 'lagrangianDataGlobal.m' or 'lagrangianDataPlanar.m' Without Further Processing
% ----
% Usage: [deltaT, particleData, particleProps] = storeLagrangianData(caseFolder, format)
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Desired Lagarangian Data Type Stored as String
%                        'global' or 'planar'


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Global and PODprobe Directory Identification
% v1.2 - Added Support for 'lagrangianDataPlanar.m'


%% Main Function

function [particleData, particleProps] = storeLagrangianData(caseFolder, format)

    clc;
    
    [timeDirs, ~] = timeDirectories(caseFolder, 'global');

    disp(' ');

    if strcmp(format, 'global')
        [particleData, particleProps] = lagrangianDataGlobal(caseFolder, timeDirs);
    elseif strcmp(format, 'planar')
        [particleData, particleProps] = lagrangianDataPlanar(caseFolder, timeDirs);
    else
        error('Invalid Lagrangian Data Format')
    end

end
