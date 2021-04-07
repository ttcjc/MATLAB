%% Lagrangian Data Storage v1.0
% ----
% Executes 'lagrangianData.m' Without Further Processing
% ----
% Usage: [particleData, particleProps] = storeLagrangianData(caseFolder)
%        'caseFolder' -> Case Path Stored as String


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [particleData, particleProps] = storeLagrangianData(caseFolder)

    [timeDirs, ~] = timeDirectories(caseFolder);

    disp(' ');
    disp(' ');

    [particleData, particleProps] = lagrangianData(caseFolder, timeDirs);

end
