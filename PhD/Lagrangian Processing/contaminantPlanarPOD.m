%% Lagrangian Contaminant Planar POD Calculator v2.0

clear variables;
close all;
clc;
% evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

% nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('=================================================');
disp ('Lagrangian Contaminant Planar POD Calculator v2.0');
disp ('=================================================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialise Case

[caseFolder, caseName, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


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
clear valid;

disp(' ');
disp(' ');


%% Acquire Contaminant Map Data

disp('Contaminant Map Data Acquisition');
disp('---------------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile(['/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/*.*'], ...
                                      'Select Map Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, ['contaminantMap/Base/', caseName])
                disp(['    Loading ''', fileName, '''...']);
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                disp('        Success');
                valid = true;
            else
                disp('    WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, ['contaminantMap/X_*/', caseName])
                disp(['    Loading ''', fileName, '''...']);
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                disp('        Success');
                valid = true;
            else
                disp('    WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
end
clear valid;

disp(' ');
disp(' ');


%% Select Mapping Options

disp('POD Options');
disp('------------');

% Select Variable of Interest
PODvar = fieldnames(mapData.mean);
PODvar = PODvar(1:(end - 1));

valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'multiple', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', PODvar);

    if ~valid
        disp(' ');
        disp('WARNING: No Mapping Variable Selected');
    end
end
clear valid;

PODvar = PODvar(index);

% Specify Map Boundaries




disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)

disp('Planar POD');
disp('-----------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;
evalc('parpool(nProc);');

disp(' ');

disp('    Initialising...');