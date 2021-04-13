%% Lagrangian Data Storage v1.1
% ----
% Executes 'lagrangianData.m' Without Further Processing
% ----
% Usage: [deltaT, particleData, particleProps] = storeLagrangianData(caseFolder)
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Required Time Directory Type Stored as String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Global and PODprobe Directory Identification


%% Main Function

function [deltaT, particleData, particleProps] = storeLagrangianData(caseFolder, format)

    clc;
    
    [timeDirs, deltaT] = timeDirectories(caseFolder, format);

    disp(' ');
    disp(' ');

    [particleData, particleProps] = lagrangianData(caseFolder, timeDirs);

end
