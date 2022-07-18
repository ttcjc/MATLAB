%% Lagrangian Particle Tracker v3.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

normalise = true; % Normalisation of Dimensions

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('=====================');
disp('Particle Tracker v3.0');
disp('=====================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Support Tracking of Specific Times Instances
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialise Case

[caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Rear-End Contamination');
disp('    B: Downstream Contaminant Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Region of Interest [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

switch format

    case 'A'
        [dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
                                                                              cloudName, false, true, ...
                                                                              true, timeDirs, deltaT, ...
                                                                              timePrecision, nProc);

    case 'B'
        [dataID, LagProps, ~, ~, LagData, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseName, ...
                                                                              cloudName, true, false, ...
                                                                              true, timeDirs, deltaT, ...
                                                                              timePrecision, nProc);
end

if normalise
    dataID = [dataID, '_Norm'];
end

disp(' ');
disp(' ');


%% Select Tracking Options

disp('Tracking Options');
disp('-----------------');

dLims = zeros(2,1);

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = [1; 120];
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        dLims(1) = inputD('Min');

        if dLims(1) == -1
            continue
        end

        dLims(2) = inputD('Max');

        if dLims(2) == -1
            continue
        end
        
        dLims = sort(dLims);
        dLims(1) = floor(dLims(1));
        dLims(2) = ceil(dLims(2));
        
        if (dLims(2) < 1) || (dLims(1) > 120)
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            continue
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if normalise
    dataID = insertBefore(dataID, '_Norm', ['_D', num2str(dLims(1)), '_D', num2str(dLims(2))]);
else
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Enable Particle Position Interpolation? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        interp = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        interp = true;
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');
























%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        D = -1;
    end
    
end