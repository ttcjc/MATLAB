%% Probe Data Initialisation v1.1
% ----
% Initialisation of OpenFOAM v7 Probe Data for Further Processing
% ----
% Usage: [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialiseProbeData(field, normalise, nProc);
%        'field'     -> Desired Field Stored as String
%        'normalise' -> Normalise Geometry Dimensions [True/False]
%        'nProc'      -> Number of Processors Used for Parallel Collation


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Update to Support Additional Versatility of 'velocityProcessing.m'


%% Supported OpenFOAM Cases

% Run_Test
% Windsor_2022


%% Supported Fields

% Pressure: 'p'
% Velocity: 'U'


%% Main Function

function [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialiseProbeData(field, normalise, nProc)
    
    % Select Case
    disp('Case Selection');
    disp('---------------');
    
    caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
    namePos = max(strfind(caseFolder, '/')) + 1;
    
    disp(' ');
    
    disp(['Case: ', caseFolder]);
    
    % Confirm Support
    if ~contains(caseFolder, ["Run_Test", "Windsor"])
        error('Invalid Case Directory (Unsupported Case Type)');
    end

    disp(' ');
    
    % Confirm Data Availability and Identify Time Directories
    if strcmp(field, 'p')
        probeType = 'probesPressure';
        [timeDirs, ~] = timeDirectories(caseFolder, probeType);
    elseif strcmp(field, 'U')
        probeType = 'probesVelocity';
        [timeDirs, ~] = timeDirectories(caseFolder, probeType);
    else
        error('Invalid Field (Available Options: ''p'' or ''U'')');
    end
    
    disp(' ');
    disp(' ');

    % Acquire Probe Data
    disp('Probe Data Acquisition');
    disp('-----------------------');
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Load Saved Probe Data? [y/n]: ', 's');
    
        if selection == 'n' | selection == 'N' %#ok<OR2>
            disp(' ');
            data = readProbeData(caseFolder, timeDirs, field, nProc);
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2> 
            [fileName, filePath] = uigetfile(['~/Data/Numerical/MATLAB/probeData/', caseFolder(namePos(end):end), '/*.*'], ...
                                             'Select Probe Data');
    
            if contains(filePath, [caseFolder(namePos(end):end), '/', probeType])
                disp(['    Loading: ', fileName]);
                load([filePath, fileName], 'data');
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

    disp(' ');
    disp(' ');
    
    % Select Relevant Geometry and Define Bounding Box
    [geometry, xDims, yDims, zDims, precision] = selectGeometry(normalise);

    % Normalise Data Dimensions
    if normalise
        
        if contains(caseFolder, ["Run_Test", "Windsor"])
            data.position(:,1) = round(data.position(:,1) / 1.044, precision);
            data.position(:,2) = round(data.position(:,2) / 1.044, precision);
            data.position(:,3) = round(data.position(:,3) / 1.044, precision);
        end
        
    end

end