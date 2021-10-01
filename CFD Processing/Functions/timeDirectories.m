%% Time Directory Reader v1.1
% ----
% Collates OpenFOAM v7 Time Directories
% ----
% Usage: [timeDirs, deltaT] = timeDirectories(caseFolder, format)
%        'caseFolder' -> Case Path Stored as String
%        'format'     -> Required Time Directory Type Stored as String
%                        'global' or 'PODprobe'

%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Support for Global and PODprobe Directory Identification


%% Main Function

function [timeDirs, deltaT] = timeDirectories(caseFolder, format)

    % Confirm Valid Case
    disp('Analysing Case Structure:');
    
    if contains(format, 'global')
        
        if exist([caseFolder, '/0'], 'dir')
            timeDirs = dir(caseFolder);
        else
            error('Invalid Case Directory (Unexpected File Structure)');
        end
        
    elseif contains(format, 'PODprobe')
        
        if exist([caseFolder, '/postProcessing/probesPOD'], 'dir')
            timeDirs = dir([caseFolder, '/postProcessing/probesPOD']);
        else
            error('Invalid Case Directory (No Probe Data Found)');
        end
        
    else
        error('Invalid Format (Available Options: ''global'' or ''PODprobe''')
    end

    % Identify Time Directories
    i = 1;
    while i <= height(timeDirs)

        if str2double(timeDirs(i,1).name) == 0 || isnan(str2double(timeDirs(i,1).name))
            timeDirs(i,:) = [];
        else
            i = i + 1;
        end

    end

    deltaT = str2double(timeDirs(end,1).name) - str2double(timeDirs(end-1,1).name);
    deltaT = str2double(num2str(deltaT));
    
    disp(['    Identified ', num2str(size(timeDirs,1)), ' Time Directories']);
    disp(['    ', char(916), 'T = ' num2str(deltaT), 's']);
    
end