%% Lagrangian Data Load v1.2
% ----
% Load Lagrangian Data Acquired Using OpenFOAM v7


%% Preamble

run preamble;

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

disp('=========================');
disp('Lagrangian Data Load v1.2');
disp('=========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Minor Formatting Updates
% v1.2 - Minor Update to Shift Preamble Into Separate Script


%% Select Data

% Select Case
disp('Case Selection');
disp('---------------');

caseFolder = uigetdir('~/OpenFOAM/', 'Select Case');

campaignID = caseFolder((strfind(caseFolder, 'results/') + 8):(max(strfind(caseFolder, '/')) - 1));
caseID = caseFolder((max(strfind(caseFolder, '/')) + 1):end);

disp(' ');

disp(['Case: ', caseID]);

% Confirm Support
if ~contains(caseID, ["Run_Test", "Windsor"])
    error('Invalid Case Directory (Unsupported Case Type)');
end

disp(' ');

% Confirm Basic Data Availability and Identify Time Directories
[timeDirs, deltaT, timePrecision] = timeDirectories(caseFolder, 'global');

% Select Data Types
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

if surface || plane || volume
    maxNumCompThreads(nProc);
    
    [dataID, LagProps, LagDataSurface, LagDataPlane, ...
     LagDataVolume, sampleInt, dataFormat] = initialiseLagData(saveLoc, caseFolder, campaignID, ...
                                                               caseID, cloudName, surface, plane, ...
                                                               volume, timeDirs, deltaT, ...
                                                               timePrecision, nProc);
end
