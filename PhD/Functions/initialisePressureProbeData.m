%% Pressure Probe Data Initialisation v1.0
% ----
% Initialisation of OpenFOAM v7 Surface Pressure Probe Data for Further Processing
% ----
% Usage: [caseFolder, probeData, timePrecision, geometry, ...
%         xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData(normalise, nProc);
%        'normalise' -> Normalise Dimensions [True/False]
%        'nProc'     -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Main Function

function [caseName, probeData, timePrecision, geometry, ...
          xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(normalise, nProc)
    
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
    [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'probesPressure');
    
    disp(' ');
    disp(' ');

    % Acquire Volumetric Probe Data
    disp('Surface Probe Data Acquisition');
    disp('-------------------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Saved Probe Data? [y/n]: ', 's');
    
        if selection == 'n' | selection == 'N' %#ok<OR2>
            disp(' ');
            probeData = readProbeData(caseFolder, caseName, timeDirs, deltaT, timePrecision, 'probesPressure', nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/probesPressure/*.mat'], ...
                                             'Select Probe Data');
    
            if contains(filePath, [caseName, '/probesPressure'])
                disp(['    Loading ''', fileName, '''...']);
                probeData = load([filePath, fileName], 'probeData').probeData;
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
        
        if contains(caseName, ["Run_Test", "Windsor"])
            probeData.position = round((probeData.position / 1.044), spacePrecision);
            
            probeData.planePosition = round((probeData.planePosition / 1.044), spacePrecision);
        end
        
    end
    
    probeData = orderfields(probeData);
    
end