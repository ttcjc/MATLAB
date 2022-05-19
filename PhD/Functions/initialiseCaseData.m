%% OpenFOAM Case Initialisation v2.0
% ----
% Initialisation of OpenFOAM v7 Case Data for Further Processing
% ----
% Usage: [caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise);
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Balance and Upstream Windsor Case Variants
% v1.2 - Added Support for Global and PODprobe Directory Identification
% v2.0 - Substantial Rewrite to Accommodate New Data Formats


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Main Function

function [caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseName = caseFolder(namePos:end);

    disp(' ');

    disp(['Case: ', caseFolder]);

    % Confirm Support
    if ~contains(caseFolder, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    disp(' ');
    
    % Confirm Basic Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');
    
    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);
    
end
    