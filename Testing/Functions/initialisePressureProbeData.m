%% Pressure Probe Data Initialisation v1.1
% ----
% Initialisation of OpenFOAM v7 Surface Pressure Probe Data for Further Processing
% ----
% Usage: [caseFolder, probeData, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(saveLocation, normalise, nProc);
%
%        'saveLocation'  -> Start of File Path, Stored as a String
%        'normalise' -> Normalise Dimensions [True/False]
%        'nProc'     -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update To Improve Consistency of Structures Across Repository


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022
% Windsor_Upstream_2023
% Windsor_fullscale


%% Main Function

function [caseID, dataID, probeData, sampleInterval, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(saveLocation, normalise, nProc)
    
    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseID = caseFolder(namePos:end);
    
    disp(' ');
    
    disp(['Case: ', caseID]);
    
    % Confirm Support
    if ~contains(caseID, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    disp(' ');
    
    % Confirm Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'probesPressure');
    
    disp(' ');
    disp(' ');

    % Acquire Surface Probe Data
    disp('Surface Probe Data Acquisition');
    disp('-------------------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Saved Probe Data? [y/n]: ', 's');
    
        if selection == 'n' | selection == 'N' %#ok<OR2>
            disp(' ');
            [dataID, probeData, sampleInterval] = readProbeData(saveLocation, caseFolder, caseID, timeDirs, deltaT, ...
                                                                timePrecision, 'probesPressure', nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/probeData/', caseID, '/probesPressure/*.mat'], ...
                                             'Select Probe Data');
    
            if contains(filePath, [caseID, '/probesPressure'])
                disp(['    Loading ''', fileName, '''...']);
                
                dataID = load([filePath, fileName], 'dataID').dataID;
                probeData = load([filePath, fileName], 'probeData').probeData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                
                disp('        Success');
                
                valid = true;
            else
                disp('    WARNING: Invalid File Selection');
                clear fileName filePath;
            end
    
        else
            disp('    WARNING: Invalid Entry');
            clear fileName filePath;
        end
    
    end
    clear valid;
    
    % Populate Planar Information
    probeData.planeOrientation = 'YZ';
    probeData.planePosition = probeData.position(1,1);
    
    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if contains(caseID, ["Run_Test", "Windsor"])
            probeData.position = round((probeData.position / 1.044), spacePrecision);
            
            probeData.planePosition = round((probeData.planePosition / 1.044), spacePrecision);
        end
        
    end
    
    probeData = orderfields(probeData);
    
end