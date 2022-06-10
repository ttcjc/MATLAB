%% POD Temporal Analysis Tool v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===============================');
disp ('POD Temporal Analysis Tool v1.0');
disp ('===============================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Analysis Type

disp('Analysis Type');
disp('--------------');

disp(' ');

disp('Possible Analysis Types:');
disp('    A: Single-Field Fourier Transform');
disp('    B: Dual-Field Temporal Comparison');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Analysis Type [A/B]: ', 's');

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


%% Load POD Data

disp('POD Data Acquisition');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/*.mat', ...
                                      'Select POD Data');

    if contains(filePath, 'POD')
        disp(['Loading ''', fileName, '''...']);
        temporalData.A.dataID = load([filePath, fileName], 'dataID').dataID;
        temporalData.A.PODdata = load([filePath, fileName], 'PODdata').PODdata;
        disp('    Success');
        
        valid = true;
    else
        disp('WARNING: Invalid File Selection');
        clear fileName filePath;
    end
    
end
clear valid;

if contains(filePath, 'planarContaminantPOD')
    namePos = strfind(filePath, '/');
    fieldName = filePath((namePos(end - 1) + 1):(namePos(end) - 1));
elseif contains(filePath, 'planarPressurePOD')
    fieldName = 'pressure';
elseif contains(filePath, 'planarVelocityPOD')
    fieldName = 'velocity';
else
    error('Unknown POD Format');
end

temporalData = renameStructField(temporalData, 'A', fieldName);

switch format
    
    case 'B'
        
        valid = false;
        while ~valid
            disp(' ');
            [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/*.mat', ...
                                             'Select POD Data');
                                         
            if contains(filePath, 'POD')
                disp(['Loading ''', fileName, '''...']);
                temporalData.B.dataID = load([filePath, fileName], 'dataID').dataID;
                temporalData.B.PODdata = load([filePath, fileName], 'PODdata').PODdata;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        end
        clear valid;
        
        if contains(filePath, 'planarContaminantPOD')
            namePos = strfind(filePath, '/');
            fieldName = filePath((namePos(end - 1) + 1):(namePos(end) - 1));
        elseif contains(filePath, 'planarPressurePOD')
            fieldName = 'pressure';
        elseif contains(filePath, 'planarVelocityPOD')
            fieldName = 'velocity';
        else
            error('Unknown POD Format');
        end
        
        temporalData = renameStructField(temporalData, 'B', fieldName);

end

disp(' ');
disp(' ');


%% Perform Fourier Transfom

fields = fieldnames(temporalData);

for i = 1:height(fields)
end