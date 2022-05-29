%% Probe Data Initialisation v1.2
% ----
% Initialisation of OpenFOAM v7 Probe Data for Further Processing
% ----
% Usage: [caseFolder, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseProbeData(field, planar, normalise, nProc);
%        'field'     -> Desired Field Stored as String
%        'planar'    -> Extract Planar Data [True/False]
%        'normalise' -> Normalise Dimensions [True/False]
%        'nProc'     -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'
% v1.2 - Added Support for Volumetric or Planar Data Extraction


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseProbeData(field, planar, normalise, nProc)
    
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
    if strcmp(field, 'p')
        probeType = 'probesPressure';
        [timeDirs, ~, timePrecision] = timeDirectories(caseFolder, probeType);
    elseif strcmp(field, 'U')
        probeType = 'probesVelocity';
        [timeDirs, ~, timePrecision] = timeDirectories(caseFolder, probeType);
    else
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end
    
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
            data = readProbeData(caseFolder, timeDirs, timePrecision, field, nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/probeData/', caseName, '/*.*'], ...
                                             'Select Probe Data');
    
            if contains(filePath, [caseName, '/', probeType])
                disp(['    Loading ''', fileName, '''...']);
                data = load([filePath, fileName], 'data').data;
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
    if planar
        data = extractPlanarProbeData(caseName, data);

        disp(' ');
        disp(' ');
    end
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if planar
            
            if contains(caseName, ["Run_Test", "Windsor"])
                planes = fieldnames(data);
                
                for i = 1:height(planes)
                    data.(planes{i}).position = round((data.(planes{i}).position / 1.044), spacePrecision);
                    
                    data.(planes{i}).xLims = round((data.(planes{i}).xLims / 1.044), spacePrecision);
                    data.(planes{i}).yLims = round((data.(planes{i}).yLims / 1.044), spacePrecision);
                    data.(planes{i}).zLims = round((data.(planes{i}).zLims / 1.044), spacePrecision);
                end
                
            end
            
        else
            
            if contains(caseName, ["Run_Test", "Windsor"])
                data.position(:,1) = round((data.position / 1.044), spacePrecision);
            end
            
        end
        
    end
    
end