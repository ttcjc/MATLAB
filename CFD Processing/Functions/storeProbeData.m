%% POD Probe Data Storage v1.0
% ----
% Executes 'PODprobeData.m' Without Further Processing
% ----
% Usage: [deltaT, probeData] = storePODprobeData(caseFolder, format)
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Required Time Directory Type Stored as String
%                        'PODprobe'


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [probeData] = storePODprobeData(caseFolder, format)

    clc;

    [timeDirs, ~] = timeDirectories(caseFolder, format);

    disp(' ');
    disp(' ');

    [probeData] = PODprobeData(caseFolder, timeDirs);

end