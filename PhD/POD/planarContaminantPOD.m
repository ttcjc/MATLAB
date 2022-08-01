%% Planar Lagrangian Contaminant POD Calculator v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('======================================');
disp('Planar Contaminant POD Calculator v2.0');
disp('======================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base Contamination');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B]: ', 's');

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


%% Acquire Contaminant Map

disp('Contaminant Map Acquisition');
disp('----------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/contaminantMap/*.mat'], ...
                                      'Select Map Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/base')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/X_')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
end
clear valid;

namePos = strfind(filePath, '/');
caseName = filePath((namePos(end - 2) + 1):(namePos(end - 1) - 1));

timePrecision = strfind(fileName, '_T') - 3;

switch format
    
    case 'B'
        planePos = filePath((namePos(end - 1) + 1):(namePos(end) - 1));

end

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Select POD Options

disp('POD Options');
disp('------------');

% Select Variable of Interest
PODvar = fieldnames(mapData.mean);
PODvar = PODvar(1:(end - 1));

valid = false;
while ~valid
    disp(' ');
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', PODvar);

    if ~valid
        disp('WARNING: No Mapping Variable Selected');
    end
end
clear valid;

PODvar = PODvar{index};

disp(['Variable of Interest: ', PODvar]);

disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)

disp('Planar Proper Orthogonal Decomposition');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Specify Map Boundaries
switch format
    
    case 'A'
        % Identify Model Base
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            
            if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
                break;
            end
            
            if i == height(parts)
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end
            
        end
    
        geoPoints = geometry.(parts{i}).vertices;
        basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.0025, 'jointType', 'square');
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        clear basePoints basePoly;
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
    
        if contains(caseName, ["Run_Test", "Windsor"])
            mapPerim = [];
            
            xLimsData = mapData.positionGrid(1,1);
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.739];
            
            if normalise
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end
            
        end
        
end

% Initialise POD Variables
PODdata.positionGrid = mapData.positionGrid;
PODdata.time = mapData.inst.time;
PODdata.(PODvar).mean = mapData.mean.(PODvar);
PODdata.(PODvar).inst = mapData.inst.(PODvar);

clear mapData;

disp(' ');

disp('    Calculating Instantaneous Field Fluctuations...');

% Initialise Progress Bar
wB = waitbar(0, ['Calculating Instantaneous ''', PODvar, ''' Fluctuations'], 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Instantaneous Variable Fluctuations
PODdata.(PODvar).prime = cell(height(PODdata.time),1);

for i = 1:height(PODdata.time)
    PODdata.(PODvar).prime{i} = PODdata.(PODvar).inst{i} - PODdata.(PODvar).mean;
    
    waitbar((i / height(PODdata.time)), wB);
end

delete(wB);

disp(' ');

% Perform Planar Snapshot POD
switch format

    case 'A'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, 'scalar', 'Base');

    case 'B'
        [fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, PODvar, 'scalar', planePos);

end

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Mode Presentation Options

disp('Mode Presentation Options');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Individual Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotModes = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = true;
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue;
        else
            nModes = sort(nModes);
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Present POD Modes

disp('Mode Presentation');
disp('------------------');

disp(' ');

if plotModes
    % Define Plot Limits
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.52725];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = PODdata.positionGrid;
    cMap = cool2warm(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-1; 1];

    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.phi_mode(:,i), -1, 1);

        switch format

            case 'A'
                figName = ['Base_POD_', PODvar, '_M', num2str(i)];

            case 'B'
                figName = [planePos, '_POD_', PODvar, '_M', num2str(i)];

        end

        CoM = [];
        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '\it{%}'];

        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    end
    
else
    disp('    Skipping Mode Presentation');
end

disp(' ');
disp(' ');


%% Save POD Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'PODdata', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'PODdata', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Select Reconstruction Options

disp('Reconstruction Options');
disp('-----------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Perform Field Reconstruction Using N Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        return;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue;
        else
            nModes = sort(nModes);
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Perform Field Reconstruction

disp('Field Reconstruction');
disp('---------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Initialise Reconstruction Variables
reconData.positionGrid = PODdata.positionGrid;
reconData.time = PODdata.time;
reconData.(PODvar).mean = PODdata.(PODvar).mean;
reconData.(PODvar).inst = cell(Nt,1);

for i = 1:Nt
    reconData.(PODvar).inst{i} = reconData.(PODvar).mean;
end

disp(' ');

% Perform Field Reconstruction
reconData = reconstructPOD(reconData, PODdata, PODvar, nModes, Ns, Nt, 'scalar', true);

if any(strcmp(PODvar, {'mass', 'massNorm'}))
    disp(' ');

    disp('    Calculating Reconstructed Centre of Mass...');

    % Initialise Progress Bar
    wB = waitbar(0, 'Calculating Reconstructed Centre of Mass', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Calculate Reconstructed CoM
    reconData.(PODvar).CoM = cell(Nt,1);
    
    for i = 1:Nt
        reconData.(PODvar).CoM{i} = zeros(1,3);
    
        reconData.(PODvar).CoM{i}(1) = reconData.positionGrid(1,1);
        reconData.(PODvar).CoM{i}(2) = sum(reconData.(PODvar).inst{i} .* reconData.positionGrid(:,2)) / sum(reconData.(PODvar).inst{i});
        reconData.(PODvar).CoM{i}(3) = sum(reconData.(PODvar).inst{i} .* reconData.positionGrid(:,3)) / sum(reconData.(PODvar).inst{i});
        
        waitbar((i / Nt), wB);
    end
    
    delete(wB);
end

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Reconstruction Presentation Options

disp('Reconstruction Presentation Options');
disp('------------------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Reconstructed Field? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotRecon = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotRecon = true;
        nFrames = inputFrames(Nt);
        
        if nFrames == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


disp(' ');
disp(' ');


%% Present Reconstruction

disp('Reconstruction Presentation');
disp('----------------------------');

disp(' ');

if plotRecon
    disp('    Presenting Reconstructed Field...');

    % Define Plot Limits
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.52725];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = reconData.positionGrid;
    cMap = viridis(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    
    if any(strcmp(PODvar, {'d10', 'd20', 'd30', 'd32'}))
        cLims = dLims;
    elseif strcmp(PODvar, 'massNorm')
        cLims = [0; 1];
    else
        cLims = [0; max(cellfun(@max, reconData.(PODvar).inst))];
    end

    figHold = fig;

    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.(PODvar).inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        
        switch format
            
            case 'A'
                figName = ['Base_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
            
            case 'B'
                figName = [planePos, '_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
        
        end
        
        if any(strcmp(PODvar, {'mass', 'massNorm'}))
            CoM = reconData.(PODvar).CoM{i};
        else
            CoM = [];
        end
        
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    end
    
else
    disp('    Skipping Reconstruction Presentation');
end

disp(' ');
disp(' ');


%% Save Reconstruction Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'dataID', 'reconData', 'nModes', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar, '/', dataID, '_', mat2str(nModes), '.mat'], ...
                     'dataID', 'reconData', 'nModes', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function modes = inputModes(Nt)

    modes = str2num(input('    Input Desired Modes [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes > Nt)
        disp('        WARNING: Invalid Entry');
        modes = -1;
    end

end


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end
