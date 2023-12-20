%% OpenFOAM Case Initialisation v2.3
% ----
% Initialisation of OpenFOAM v7 Case Data for Further Processing
% ----
% Usage: [caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData;


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Balance and Upstream Windsor Case Variants
% v1.2 - Added Support for Global and PODprobe Directory Identification
% v2.0 - Substantial Rewrite to Accommodate New Data Formats
% v2.1 - Added Support for Full-Scale Windsor Model Simulations
% v2.2 - Added 'campaignID' Output
% v2.3 - Removed Geometry Normalisation (Should Be Done After Processing)


%% Supported OpenFOAM Campaigns

% Windsor_Upstream_2023
% Windsor_fullScale


%% Main Function

function [caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData

    % Select Case
    disp('Case Selection');
    disp('---------------');

    caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');
    
    campaignID = caseFolder((strfind(caseFolder, 'results/') + 8):(max(strfind(caseFolder, '/')) - 1));
    caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

    disp(' ');

    disp(['Case: ', campaignID, ', ', caseID]);

    % Confirm Support
    if ~strcmp(campaignID, 'Windsor_Upstream_2023') && ~strcmp(campaignID, 'Windsor_fullScale')
        error('Invalid Case Directory (Unsupported Case Type)');
    end
    
    disp(' ');
    
    % Confirm Basic Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');
    
    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry;
    
end
