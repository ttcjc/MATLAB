%% OpenFOAM Case Initialisation v2.1
% ----
% Initialisation of OpenFOAM v7 Case Data for Further Processing
% ----
% Usage: [caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision, normalise, normLength] = initialiseCaseData(normalise);
%        'normalise' -> Normalise Dimensions [True/False]


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Balance and Upstream Windsor Case Variants
% v1.2 - Added Support for Global and PODprobe Directory Identification
% v2.0 - Substantial Rewrite to Accommodate New Data Formats
% v2.1 - Added Support for Full-Scale Windsor Model Simulations


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Main Function

function [caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision, normalise, normLength] = initialiseCaseData(normalise)

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseID = caseFolder(namePos:end);

    disp(' ');

    disp(['Case: ', caseID]);

    % Confirm Support
    if ~contains(caseID, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    disp(' ');
    
    % Confirm Basic Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');
    
    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise, normLength] = selectGeometry(normalise);
    
end