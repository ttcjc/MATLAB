%% POD Probe Data Storage v1.0
% ----
% Executes 'PODprobeData.m' Without Further Processing
% ----
% Usage: [probeData] = storePODprobeData(caseFolder)
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Required Time Directory Type Stored as String


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [deltaT, probeData] = storePODprobeData(caseFolder, format)

    clc;

    [timeDirs, deltaT] = timeDirectories(caseFolder, format);

    disp(' ');
    disp(' ');

    [probeData] = PODprobeData(caseFolder, timeDirs);

end
