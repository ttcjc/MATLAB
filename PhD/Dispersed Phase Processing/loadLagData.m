%% Lagrangian Data Load v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

% saveLocation = '/mnt/Processing/Data';
saveLocation = '~/Data';

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

disp('=========================');
disp('Lagrangian Data Load v1.0');
disp('=========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Data

% Select Case
disp('Case Selection');
disp('---------------');

caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');

namePos = max(strfind(caseFolder, '/')) + 1;
caseName = caseFolder(namePos:end);

disp(' ');

disp(['Case: ', caseName]);

% Confirm Support
if ~contains(caseName, ["Run_Test", "Windsor"])
    error('Invalid Case Directory (Unsupported Case Type)');
end

disp(' ');

% Confirm Basic Data Availability and Identify Time Directories
[timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');

% Select Data Types
valid = false;
while ~valid
    disp(' ');
    selection = input('Collate Plane Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plane = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plane = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    selection = input('Collate Surface Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        surface = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        surface = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    selection = input('Collate Volume Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        volume = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        volume = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Collate Data

[dataID, LagProps, LagDataPlane, LagDataSurface, ...
     LagDataVolume, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, cloudName, ...
                                                        plane, surface, volume, ...
                                                        timeDirs, deltaT, timePrecision, nProc);