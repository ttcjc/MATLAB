%% Time Directory Reader v1.2
% ----
% Collates OpenFOAM v7 Time Directories
% ----
% Usage: [timeDirs, deltaT] = timeDirectories(caseFolder, format);
%
%        'caseFolder' -> Case Path, Stored as a String
%        'format'     -> Required Time Directory Type, Stored as a String


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for 'global' and 'probe' Directory Identification
% v1.2 - Minor Update to Sort Time Directories (Allowing for 'timeDirs' > 10 s)


%% Supported Formats

% Primary Time Directories: 'global'
% Pressure Probe Time Directories: 'pProbes'
% Velocity Probe Time Directories: 'uProbes'


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

        case 'pProbes'
            
            if exist([caseFolder, '/postProcessing/probesPressure'], 'dir')
                timeDirs = dir([caseFolder, '/postProcessing/probesPressure']);
            else
                error('Invalid Case Directory (No Pressure Probe Data Found)');
            end
        
        case 'uProbes'
            
            if exist([caseFolder, '/postProcessing/probesVelocity'], 'dir')
                timeDirs = dir([caseFolder, '/postProcessing/probesVelocity']);
            else
                error('Invalid Case Directory (No Velocity Probe Data Found)');
            end
        
        otherwise
            error('Invalid Format (Available Options: ''global'', ''pProbes'' or ''uProbes'')');

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
    
    % Sort Time Directories
    [~, index] = sort(str2double({timeDirs.name}));
    timeDirs = timeDirs(index);

    % Determine Delta T
    deltaT = str2double(timeDirs(end).name) - str2double(timeDirs(end - 1).name);
    timePrecision = width(extractAfter(num2str(deltaT, 8), '.'));
    deltaT = round(deltaT, timePrecision);
    
    disp(['        Identified ', num2str(height(timeDirs)), ' Time Directories']);
    disp(['        ', char(916), 'T = ' num2str(deltaT, ['%.', num2str(timePrecision), 'f']), 's']);
    
end
