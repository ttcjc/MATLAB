%% Velocity Probe Data Initialisation v1.2
% ----
% Initialisation of OpenFOAM v7 Volumetric Velocity Probe Data for Further Processing
% ----
% Usage: [caseFolder, probeData, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(planar, normalise, nProc);
%        'format'    -> Data Processing Format Stored as String
%        'normalise' -> Normalise Dimensions [True/False]
%        'nProc'     -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v1.2 - Added Support for Volumetric or Planar Data Extraction


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Processing Formats

% Volumetric: 'volume'
% Planar:     'planar'
% Planar POD: 'planarPOD'


%% Main Function

function [caseName, dataID, probeData, sampleInterval, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(format, normalise, nProc)
    
    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
    
    namePos = max(strfind(caseFolder, '/')) + 1;
    caseName = caseFolder(namePos:end);
    
    disp(' ');
    
    disp(['Case: ', caseName]);
    
    % Confirm Support
    if ~contains(caseName, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    disp(' ');
    
    % Confirm Data Availability and Identify Time Directories
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'probesVelocity');
    
    disp(' ');
    disp(' ');

    % Acquire Volumetric Probe Data
    disp('Volumetric Probe Data Acquisition');
    disp('----------------------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Saved Probe Data? [y/n]: ', 's');
    
        if selection == 'n' | selection == 'N' %#ok<OR2>
            disp(' ');
            [dataID, probeData, sampleInterval] = readProbeData(caseFolder, caseName, timeDirs, deltaT, ...
                                                                timePrecision, 'probesVelocity', nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/probesVelocity/*.mat'], ...
                                             'Select Probe Data');
    
            if contains(filePath, [caseName, '/probesVelocity'])
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
    
    disp(' ');
    disp(' ');
        
    % Extract Planar Probe Data
    if strcmp(format, 'planar')
        probeData = extractPlanarVelocityProbeData(caseName, probeData, true);
    elseif strcmp(format, 'planarPOD')
        probeData = extractPlanarVelocityProbeData(caseName, probeData, false);
    end
    
    if strcmp(format, 'planar') || strcmp(format, 'planarPOD')
        disp(' ');
        disp(' ');
    end
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if strcmp(format, 'planar') || strcmp(format, 'planarPOD')
            
            if contains(caseName, ["Run_Test", "Windsor"])
                planes = fieldnames(probeData);
                
                for i = 1:height(planes)
                    probeData.(planes{i}).position = round((probeData.(planes{i}).position / 1.044), spacePrecision);
                    
                    probeData.(planes{i}).planePosition = round((probeData.(planes{i}).planePosition / 1.044), spacePrecision);
                    probeData.(planes{i}).xLims = round((probeData.(planes{i}).xLims / 1.044), spacePrecision);
                    probeData.(planes{i}).yLims = round((probeData.(planes{i}).yLims / 1.044), spacePrecision);
                    probeData.(planes{i}).zLims = round((probeData.(planes{i}).zLims / 1.044), spacePrecision);
                end
                
            end
            
        else
            
            if contains(caseName, ["Run_Test", "Windsor"])
                probeData.position = round((probeData.position / 1.044), spacePrecision);
            end
            
        end
        
    end
    
end