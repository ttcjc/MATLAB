%% Planar Pressure Processing v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===============================');
disp ('Planar Pressure Processing v4.0');
disp ('===============================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Support ParaView, Probe and Experimental Planar Data


%% Select Data Format

disp('Data Format');
disp('------------');

disp(' ');

disp('Compatible Data Formats:');
disp('    A: ParaView Planar Data');
disp('    B: Probe Planar Data');
disp('    C: Experimental Planar Data (Varney)');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Data Format [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

disp (' ');
disp (' ');


%% Initialisation

switch format

    case 'A'
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialisePVdata('p', normalise);

    case 'B'
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseProbeData('p', false, normalise, nProc);    

    case 'C'
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseExpData('p', normalise);

end


%% Data Formatting