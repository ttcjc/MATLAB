%% Time Directory Reader v1.1
% ----
% Collates OpenFOAM v7 Time Directories
% ----
% Usage: [timeDirs, deltaT] = timeDirectories(caseFolder, format);
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Required Time Directory Type Stored as String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for 'global' and 'probe' Directory Identification


%% Supported Formats

% Primary Time Directories: 'global'
% Pressure Probe Time Directories: 'probesPressure'
% Velocity Probe Time Directories: 'probesVelocity'


%% Main Function

function [timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, format)

    % Confirm Data Availability
    disp('    Analysing Directory Structure...');

    switch format

        case 'global'
            
            if exist([caseFolder, '/0'], 'dir')
                timeDirs = dir(caseFolder);
            else
                error('Invalid Case Directory (Unexpected File Structure)');
            end

        case 'probesPressure'
            
            if exist([caseFolder, '/postProcessing/probesPressure'], 'dir')
                timeDirs = dir([caseFolder, '/postProcessing/probesPressure']);
            else
                error('Invalid Case Directory (No Pressure Probe Data Found)');
            end
        
        case 'probesVelocity'
            
            if exist([caseFolder, '/postProcessing/probesVelocity'], 'dir')
                timeDirs = dir([caseFolder, '/postProcessing/probesVelocity']);
            else
                error('Invalid Case Directory (No Velocity Probe Data Found)');
            end
        
        otherwise
            error('Invalid Format (Available Options: ''global'', ''probesPressure'' or ''probesVelocity'')');

    end

    % Identify Time Directories
    i = 1;
    while i <= height(timeDirs)

        if str2double(timeDirs(i).name) == 0 || isnan(str2double(timeDirs(i).name))
            timeDirs(i) = [];
        else
            i = i + 1;
        end

    end
    clear i;

    deltaT = str2double(timeDirs(end).name) - str2double(timeDirs(end - 1).name);
    timePrecision = width(extractAfter(num2str(deltaT, 8), '.'));
    deltaT = round(deltaT, timePrecision);
    
    disp(['        Identified ', num2str(height(timeDirs)), ' Time Directories']);
    disp(['        ', char(916), 'T = ' num2str(deltaT, ['%.', num2str(timePrecision), 'f']), 's']);
    
end
