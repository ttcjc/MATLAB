%% Lagrangian Contamination Mapper v2.0

clear variables;
close all;
clc;

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

normalise = true; % Normalisation of Dimensions

nProc = 4; % Number of Processors Used for Parallel Collation

disp ('====================================');
disp ('Lagrangian Contamination Mapper v2.0');
disp ('====================================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Updated Name and Added Time-Averaging Functionality
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats 


%% Initialise Case

[caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision] = initialiseCaseData(normalise);


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Surface Contamination (Base)');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B/C]: ', 's');

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

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

switch format
    
    case 'A'
        [LagProps, LagDataPlane, LagDataSurface, LagDataVolume] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                                                    false, true, false, ...
                                                                                    timeDirs, deltaT, timePrecision, nProc);
                                                                                
    case 'B'
        [LagProps, LagDataPlane, LagDataSurface, LagDataVolume] = initialiseLagData(caseFolder, caseName, cloudName, ...
                                                                                    true, false, false, ...
                                                                                    timeDirs, deltaT, timePrecision, nProc);

        % Select Plane of Interest
        planes = fieldnames(LagDataPlane);
        
        valid = false;
        while ~valid
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'single', ...
                                     'name', 'Select Plane of Interest', ...
                                     'listString', planes);
            
            if ~valid
                disp('WARNING: No Case Type Selected');
            end
        
        end

        LagDataPlane = LagDataPlane.(planes{index});
end

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

disp(' ');

disp('Loaded Contaminant Data Spans:');

disp(['    T = ', num2str(LagDataSurface.time(1), ['%.', num2str(timePrecision), 'f']), ' s ', ...
       '-> T = ', num2str(LagDataSurface.time(end), ['%.', num2str(timePrecision), 'f']), ' s']);
disp(['    ', char(916), 'T = ' num2str(deltaT, ['%.', num2str(timePrecision), 'f']), 's']);

valid = false;
while ~valid
    disp(' ');
    selection = input('Utilise All Available Time Instances? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        timeInsts = inputTimes(LagDataSurface.time);
        
        if timeInsts == -1
            continue
        elseif isempty(timeInsts)
            disp('        WARNING: No Lagrangian Data Available for Selected Time Instances');
            continue
        end
        
        mapData.inst.time = LagDataSurface.time(timeInsts);
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        contaminantData.time = LagDataSurface.time;

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

dLims = zeros(2,1);

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims(1) = floor(cellfun(@min, LagDataSurface.d) * 1e6); % Converts to um
        dLims(2) = ceil(cellfun(@max, LagDataSurface.d) * 1e6);
        
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
        
        if (dLims(2) < floor(cellfun(@min, LagDataSurface.d) * 1e6)) || ...
           (dLims(1) > ceil(cellfun(@max, LagDataSurface.d) * 1e6))
            disp('        WARNING: No Lagrangian Data in Selected Time Range');
            continue
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp(' ');
disp(' ');


%% Contaminant Mapping

disp('Contaminant Mapping');
disp('--------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

disp(' ');


%% Local Functions

function timeInsts = inputTimes(origTimes)

    timeInsts = str2num(input('    Input Desired Time Instances (Row Vector Form) [s]: ', 's')); %#ok<ST2NM>
    
    if any(isnan(timeInsts)) || ~isrow(timeInsts) > 1 || any(timeInsts <= 0)
        disp('        WARNING: Invalid Entry');
        timeInsts = -1;
    else
        timeInsts = find(ismember(origTimes, timeInsts));
    end

end